use sprs::CsMat;
use std::cmp::max;

mod graph;
use graph::build_graph;

mod node_math;
use node_math::{gallager_prod_exc_one, normalized_mult_exc_one,
		normalized_mult, hard_decision};

pub struct Decoder {
    pub info_pos: Vec<i32>,
    pub iter:     u32,
    pub graph:    DecoderGraph,
    pub state:    DecoderState,
    pub scratch:  DecoderScratch,
}

pub struct DecoderGraph {
    pub n:          usize,
    pub m:          usize,
    pub n_edges:    usize,
    pub cn_edges:   Vec<Vec<usize>>,
    pub vn_edges:   Vec<Vec<usize>>,
    pub cn_max_deg: usize,
    pub vn_max_deg: usize,
    pub edge_to_vn: Vec<usize>
}

pub struct DecoderState {
    pub p0_aprio:     Vec<f32>,
    pub msg_cn_to_vn: Vec<f32>,
    pub msg_vn_to_cn: Vec<f32>
}

pub struct DecoderScratch {
    pub prefix_f0: Vec<f32>,
    pub prefix_f1: Vec<f32>,
    pub suffix_f0: Vec<f32>,
    pub suffix_f1: Vec<f32>,
    pub result:    Vec<f32>
}

impl Decoder {
    // Constructor
    pub fn new(h: CsMat<u8>, info_positions: Vec<i32>) -> Self {

	let info_pos = info_positions;
	// Set default for iter, the maximum number of iterations
	let iter = 100;

	let graph   = DecoderGraph::new(h);
	let state   = DecoderState::new(&graph);
	let scratch = DecoderScratch::new(&graph);
	
	Self {
	    info_pos,
	    iter,
	    graph,
	    state,
	    scratch
	}
    }
    
    pub fn decode(&mut self, recv: &[f32], sigma: f32 ) -> Result<(), String> {
	if recv.len() != self.graph.n {
            return Err(format!(
		"decode(): recv length {} does not match code length {}",
		recv.len(), self.graph.n));
	}

	// Calculate a-priori probabilites based on channel output
	let alpha = 2.0 / (sigma * sigma);

	self.state.reset_msg();
	for i in 0..self.graph.n {
	    self.state.p0_aprio[i] = 1.0 / (1.0 + (alpha * recv[i]).exp())
	}

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
		if self.valid_cw(&vn_quant) == true {
		    println!("Valid codeword found, {} iterations used.", i);
		    break;
		}
	    }
	}
	Ok(())
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
	for edges in self.graph.cn_edges.iter() {
	    incoming.clear();
	    for &e in edges {
		incoming.push(self.state.msg_vn_to_cn[e]);
	    }
	    let deg = incoming.len();
	    let prefix_f0 = &mut self.scratch.prefix_f0[..deg];
	    let suffix_f0 = &mut self.scratch.suffix_f0[..deg];
	    let result    = &mut self.scratch.result[..deg];
	    gallager_prod_exc_one(&incoming, prefix_f0, suffix_f0, result);

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
    
    pub fn valid_cw(&self, vn_quantized: &[u8]) -> bool {
	for edges in &self.graph.cn_edges {
            let mut parity = 0u8;

            for e in edges {
		let vn = self.graph.edge_to_vn[*e];
		parity ^= vn_quantized[vn];
            }

            if parity != 0 {
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
    
	println!("Decoder properties:\nn: {}, \
		  k: {}, max iterations: {}, max dc: {}, max dv: {}",
		 self.graph.n, self.graph.n - self.graph.m, self.iter,
		 self.graph.cn_max_deg, self.graph.vn_max_deg);
	if !syst_enc {
	    println!("Using non-systematic encoding.");
	} else {
	    println!("Using systematic encoding, information positions: {:?}",
		     self.info_pos);
	}
    }
}

impl DecoderGraph {
    // Constructor
    pub fn new(h: CsMat<u8>) -> Self {

	let (m, n)  = h.shape();
	let n_edges = h.nnz();

	// Build the graph and the edges for decoding
	let (cn_edges, vn_edges, edge_to_vn) = build_graph(&h);
	let cn_max_deg = cn_edges.iter().map(|c| c.len()).max().unwrap();
	let vn_max_deg = vn_edges.iter().map(|v| v.len()).max().unwrap();

	Self {
	    n,
	    m,
	    n_edges,
	    cn_edges,
	    vn_edges,
	    cn_max_deg,
	    vn_max_deg,
	    edge_to_vn
	}
    }
}

impl DecoderState {
    // Constructor
    pub fn new(graph: &DecoderGraph) -> Self {
	let p0_aprio = vec![0.0; graph.n];
	let msg_cn_to_vn = vec![0.5; graph.n_edges]; // 0.5: first half-iter
	let msg_vn_to_cn = vec![0.0; graph.n_edges];

	Self {
	    p0_aprio,
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
    pub fn new(graph: &DecoderGraph) -> Self {
	// Allocate prefix/suffix vectors and result vector
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
