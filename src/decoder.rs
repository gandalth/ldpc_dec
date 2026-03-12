use sprs::CsMat;

mod graph;
use graph::build_graph;

mod node_math;
use node_math::{gallager_prod_exc_one, normalized_mult_exc_one,
		normalized_mult, hard_decision};

pub struct Decoder {
    pub n:            usize,
    pub k:            usize,
    pub info_pos:     Vec<i32>,
    pub iter:         u32,
    pub p0_aprio:     Vec<f32>,
    pub msg_cn_to_vn: Vec<f32>,
    pub msg_vn_to_cn: Vec<f32>,
    pub cn_edges:     Vec<Vec<usize>>,
    pub vn_edges:     Vec<Vec<usize>>,
    pub cn_max_deg:   usize,
    pub vn_max_deg:   usize,
    pub edge_to_vn:   Vec<usize>,
}

impl Decoder {
    // Constructor
    pub fn new(h: CsMat<u8>, info_positions: Vec<i32>) -> Self {

	let (k, n)  = h.shape();
	let n_edges = h.nnz();

	let p0_aprio = vec![0.0; n];
	let info_pos = info_positions;
	
	// Build the graph and the edges for decoding
	let (cn_edges, vn_edges, edge_to_vn) = build_graph(&h);
	let cn_max_deg = cn_edges.iter().map(|c| c.len()).max().unwrap();
	let vn_max_deg = vn_edges.iter().map(|v| v.len()).max().unwrap();
	
	let msg_cn_to_vn = vec![0.5; n_edges]; // Init 0.5 for first half-iter
	let msg_vn_to_cn = vec![0.0; n_edges];

	// Set default for iter, the maximum number of iterations
	let iter = 100;
	
	Self {
	    n,
	    k,
	    info_pos,
	    iter,
	    p0_aprio,
	    msg_cn_to_vn,
	    msg_vn_to_cn,
	    cn_edges,
	    vn_edges,
	    cn_max_deg,
	    vn_max_deg,
	    edge_to_vn
	}
    }

    
    pub fn decode(&mut self, recv: &[f32], sigma: f32 ) -> Result<(), String> {
	if recv.len() != self.n {
            return Err(format!(
		"decode(): recv length {} does not match code length {}",
		recv.len(), self.n));
	}

	// Calculate a-priori probabilites based on channel output
	let sigma2 = sigma * sigma;
	self.p0_aprio = recv.iter()
	    .map(|&r| 1.0 / (1.0 + (2.0 * r / sigma2).exp()))
	    .collect();
	
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
	let mut incoming = Vec::with_capacity(self.vn_max_deg);
	for (vn, edges) in self.vn_edges.iter().enumerate() {

	    incoming.clear();
	    for &e in edges {
		incoming.push(self.msg_cn_to_vn[e]);
	    }

	    let mut result = normalized_mult_exc_one(&incoming);

	    // Multiply a-prio info from channel ONCE per outgoing edge
	    let mut pair = [0.0f32; 2];
	    pair[1] = self.p0_aprio[vn];
	
	    for val in result.iter_mut() {
		pair[0] = *val;
		*val = normalized_mult(&pair);
	    }

	    for (&e, &val) in edges.iter().zip(result.iter()) {
		self.msg_vn_to_cn[e] = val;
	    }
	}
    }

    pub fn cn_update(&mut self) {
	let mut incoming = Vec::with_capacity(self.cn_max_deg);
	for edges in self.cn_edges.iter() {
	    incoming.clear();
	    for &e in edges {
		incoming.push(self.msg_vn_to_cn[e]);
	    }

	    let result = gallager_prod_exc_one(&incoming);

	    for (&e, &val) in edges.iter().zip(result.iter()) {
		self.msg_cn_to_vn[e] = val;
	    }
	}
    }

    pub fn vn_aposteriori(&self) -> Vec<f32> {
	let mut result = Vec::with_capacity(self.vn_edges.len());
	let mut incoming = Vec::with_capacity(self.vn_max_deg);
	
	for (vn, edges) in self.vn_edges.iter().enumerate() {
	    incoming.clear();
	    for &e in edges {
		incoming.push(self.msg_cn_to_vn[e]);
	    }
	    // Add apriori information to find aposteriori info per vn
	    incoming.push(self.p0_aprio[vn]);
	    
	    result.push(normalized_mult(&incoming));
	}
	return result;
    }
    
    pub fn valid_cw(&self, vn_quantized: &[u8]) -> bool {
	for edges in &self.cn_edges {
            let mut parity = 0u8;

            for e in edges {
		let vn = self.edge_to_vn[*e];
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
    
	println!("Decoder properties:\nn: {}, k: {}, max iterations: {}, max dc: {}, max dv: {}",
		 self.n, self.k, self.iter, self.cn_max_deg, self.vn_max_deg);
	if !syst_enc {
	    println!("Using non-systematic encoding.");
	} else {
	    println!("Using systematic encoding, information positions: {:?}",
		     self.info_pos);
	}
    }
}
