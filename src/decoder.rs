use std::fmt;

use std::cmp::max;

use crate::node_math::{gallager_prod_exc_one, normalized_mult_exc_one,
		       normalized_mult, hard_decision};

mod op_mode;
pub use op_mode::OpMode;

use crate::channel::Channel;
use crate::graph::Graph;

pub struct Decoder<'a> {
    pub info_pos: Vec<i32>,
    pub mode:     OpMode,
    pub iter:     u32,
    pub graph:    &'a Graph,
    pub state:    DecoderState,
    pub scratch:  DecoderScratch,
    pub result:   DecoderResult,
}

pub struct DecoderState {
    pub p0_aprio:      Vec<f32>,
    pub soft_syndrome: Vec<f32>,
    pub hard_syndrome: Vec<u8>,
    pub msg_cn_to_vn:  Vec<f32>,
    pub msg_vn_to_cn:  Vec<f32>
}

pub struct DecoderScratch {
    pub prefix_f0: Vec<f32>,
    pub prefix_f1: Vec<f32>,
    pub suffix_f0: Vec<f32>,
    pub suffix_f1: Vec<f32>,
    pub result:    Vec<f32>
}

pub struct DecoderResult {
    pub estimate: Vec<u8>,
    pub iterations: u32,
    pub success: bool
}

impl <'a> Decoder <'a> {
    // Constructor
    pub fn new(graph: &'a Graph, opmode: OpMode,
	       info_positions: Vec<i32>) -> Self {

	let info_pos = info_positions;
	let mode = opmode;
	// Set default for iter, the maximum number of iterations
	let iter = 100;

	let state   = DecoderState::new(&graph);
	let scratch = DecoderScratch::new(&graph);
	let result  = DecoderResult::new(&graph);
	
	Self {
	    info_pos,
	    mode,
	    iter,
	    graph,
	    state,
	    scratch,
	    result
	}
    }

    pub fn apply_channel<C: Channel>(&mut self, channel: &C, tx: &[C::Tx])
				     -> Result<(), String> {
            self.state.reset_msg();
            channel.apply(tx, &mut self.state)
	}

    pub fn decode(&mut self) {
	let mut i = 0u32;
	while i < self.iter {
	    self.vn_update(); // Start with vn_update to get channel info
	    self.cn_update();
	    i += 1;
	    // First full iteration before checking for valid CW.
	    // This is a small loss if a valid cw is given to the decoder.
	    if i % 5 == 1 {
		let vn_apost = self.vn_aposteriori();
		let vn_quant = hard_decision(&vn_apost);
		if self.satisfies_syndrome(&vn_quant,
					   &self.state.hard_syndrome) {
		    self.result.estimate   = vn_quant;
		    self.result.iterations = i;
		    self.result.success    = true;
		    return;
		}
	    }
	}
	// --- recompute final estimate ---
	let vn_apost = self.vn_aposteriori();
	let vn_quant = hard_decision(&vn_apost);

	self.result.estimate   = vn_quant;
	self.result.iterations = self.iter;
	self.result.success    = false;
    }
    
    pub fn vn_update(&mut self) {
	let mut incoming = Vec::with_capacity(self.graph.vn_max_deg);

	for (vn, edges) in self.graph.vn_edges.iter().enumerate() {

	    incoming.clear();
	    for &e in edges {
		incoming.push(self.state.msg_cn_to_vn[e]);
	    }

	    let deg = incoming.len();
	    let prefix_f0 = &mut self.scratch.prefix_f0[..deg];
	    let suffix_f0 = &mut self.scratch.suffix_f0[..deg];
	    let prefix_f1 = &mut self.scratch.prefix_f1[..deg];
	    let suffix_f1 = &mut self.scratch.suffix_f1[..deg];
	    let result    = &mut self.scratch.result[..deg];

	    normalized_mult_exc_one(&incoming, prefix_f0, suffix_f0,
				    prefix_f1, suffix_f1, result);
	    // Multiply a-prio info from channel ONCE per outgoing edge
	    let mut pair = [0.0f32; 2];
	    pair[1] = self.state.p0_aprio[vn];
	
	    for val in result.iter_mut() {
		pair[0] = *val;
		*val = normalized_mult(&pair);
	    }

	    for (&e, &val) in edges.iter().zip(result.iter()) {
		self.state.msg_vn_to_cn[e] = val;
	    }
	}
    }

    pub fn cn_update(&mut self) {
	let mut incoming = Vec::with_capacity(self.graph.cn_max_deg);
	for (j, edges) in self.graph.cn_edges.iter().enumerate() {
	    incoming.clear();
	    for &e in edges {
		incoming.push(self.state.msg_vn_to_cn[e]);
	    }
	    let deg = incoming.len();
	    let prefix_f0 = &mut self.scratch.prefix_f0[..deg];
	    let suffix_f0 = &mut self.scratch.suffix_f0[..deg];
	    let result    = &mut self.scratch.result[..deg];

	    let ss = self.state.soft_syndrome[j];
	    gallager_prod_exc_one(&incoming, ss, prefix_f0, suffix_f0, result);

	    for (&e, &val) in edges.iter().zip(result.iter()) {
		self.state.msg_cn_to_vn[e] = val;
	    }
	}
    }

    pub fn vn_aposteriori(&self) -> Vec<f32> {
	let mut result = Vec::with_capacity(self.graph.vn_edges.len());
	let mut incoming = Vec::with_capacity(self.graph.vn_max_deg);
	
	for (vn, edges) in self.graph.vn_edges.iter().enumerate() {
	    incoming.clear();
	    for &e in edges {
		incoming.push(self.state.msg_cn_to_vn[e]);
	    }
	    // Add apriori information to find aposteriori info per vn
	    incoming.push(self.state.p0_aprio[vn]);
	    
	    result.push(normalized_mult(&incoming));
	}
	return result;
    }
    
    pub fn satisfies_syndrome(&self, vn_quant: &[u8], syn_quant: &[u8])
			      -> bool {
	for (j, edges) in self.graph.cn_edges.iter().enumerate() {
            let mut parity = 0u8;

            for e in edges {
		let vn = self.graph.edge_to_vn[*e];
		parity ^= vn_quant[vn];
            }

            if parity !=  syn_quant[j] {
		return false;
            }
	}
	true
    }

    pub fn info(&self) {
	let syst_enc: bool; 
    
	if let Some(_max_info) = self.info_pos.iter().max() {
	    syst_enc = true;
	} else {
	    syst_enc = false;
	}
    
	println!("\nInformation:\n\
		  Decoder mode: {}, max iterations: {}\n\
		  Code properties: n: {}, k: {}, max dc: {}, max dv: {}",
		 self.mode, self.iter,
		 self.graph.n, self.graph.n - self.graph.m,
		 self.graph.cn_max_deg, self.graph.vn_max_deg);
	if !syst_enc {
	    println!("Encoding: Non-systematic.\n");
	} else {
	    println!("Encoding: Systematic, information positions: {:?}\n",
		     self.info_pos);
	}
    }
}

impl DecoderState {
    // Constructor
    pub fn new(graph: &Graph) -> Self {
	let p0_aprio      = vec![0.0; graph.n];
	let soft_syndrome = vec![0.0; graph.m];
	let hard_syndrome = vec![0;   graph.m];
	let msg_cn_to_vn  = vec![0.5; graph.n_edges]; // 0.5: first half-iter
	let msg_vn_to_cn  = vec![0.0; graph.n_edges];

	Self {
	    p0_aprio,
	    soft_syndrome,
	    hard_syndrome,
	    msg_cn_to_vn,
	    msg_vn_to_cn,
	}
    }

    pub fn reset_msg(&mut self) {
	self.msg_cn_to_vn.fill(0.5); // Init 0.5 for first half-iter
	self.msg_vn_to_cn.fill(0.0);
	// Note: filling of msg_vn_to_cn could be skipped for performance,
	// as iteration order is "vn_update first". Left in for clarity. 
    }
}


impl DecoderScratch {
    // Constructor
    pub fn new(graph: &Graph) -> Self {
	// Scratch buffers are used by kernel in node_math.rs file.
	// Resetting and filling is up to these kernels.
	let max_deg = max(graph.vn_max_deg, graph.cn_max_deg);
	let prefix_f0 = vec![1.0; max_deg];
	let prefix_f1 = vec![1.0; max_deg];
	let suffix_f0 = vec![1.0; max_deg];
	let suffix_f1 = vec![1.0; max_deg];
	let result    = vec![0.0; max_deg];

	Self {
	    prefix_f0,
	    prefix_f1,
	    suffix_f0,
	    suffix_f1,
	    result
	}
    }
}

impl DecoderResult {
    // Constructor
    pub fn new(graph: &Graph) -> Self {
	// Result of the decoding process
	let estimate   = vec![0; graph.n];
	let iterations = 0;
	let success    = false;

	Self {
	    estimate,
	    iterations,
	    success
	}
    }
}

impl fmt::Display for DecoderResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let status = if self.success { "SUCCESS" } else { "FAILURE" };

        write!(f, "{}, Iterations: {}", status, self.iterations)?;

        // Only print estimate if alternate flag is set
        if f.alternate() {
            write!(f, "Estimate: ")?;
            for bit in &self.estimate {
                write!(f, "{}", bit)?;
            }
            writeln!(f)?;
        }

        Ok(())
    }
}
