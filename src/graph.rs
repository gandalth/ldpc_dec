use sprs::CsMat;
use crate::mode::{Mode, Classic, Quantum};

use std::marker::PhantomData;

pub struct Graph<M: Mode> {
    pub n_data:         usize, // Data vn length, classic + quantum: n
    pub n_total:        usize, // Total vn length, classic: n, quantum: n+m
    pub m:              usize,
    pub n_edges:        usize,
    pub cn_edges_data:  Vec<Vec<usize>>, // Connection to data VNs
    pub cn_edges_total: Vec<Vec<usize>>, // Connection to data + extended VNs
    pub vn_edges:       Vec<Vec<usize>>,
    pub cn_max_deg:     usize,
    pub vn_max_deg:     usize,
    pub edge_to_vn:     Vec<usize>,
    _mode:              PhantomData<M>
}


impl Graph<Classic> {
    pub fn new(h: CsMat<u8>) -> Self {
        let (m, n_data) = h.shape();

        let n_total = n_data; // Classic decoding, no extra vns

	let cn_edges_data;
	let cn_edges_total;
        let vn_edges;
        let edge_to_vn;

        (cn_edges_data, cn_edges_total, vn_edges, edge_to_vn) =
                    build_classical_graph(&h);

        let cn_max_deg = cn_edges_total.iter().map(|c| c.len()).max().unwrap();
        let vn_max_deg = vn_edges.iter().map(|v| v.len()).max().unwrap();

        Self {
            n_data,
	    n_total,
            m,
            n_edges: edge_to_vn.len(),
            cn_edges_data,
	    cn_edges_total,
            vn_edges,
            cn_max_deg,
            vn_max_deg,
            edge_to_vn,
	    _mode: PhantomData,
        }
    }
}

impl Graph<Quantum> {
    pub fn new(h: CsMat<u8>) -> Self {
        let (m, n_data) = h.shape();

        let n_total = n_data + m; // Extra variable nodes for syndrome

	let cn_edges_data;
	let cn_edges_total;
        let vn_edges;
        let edge_to_vn;

        (cn_edges_data, cn_edges_total, vn_edges, edge_to_vn) =
            build_quantum_graph(&h);

        let cn_max_deg = cn_edges_total.iter().map(|c| c.len()).max().unwrap();
        let vn_max_deg = vn_edges.iter().map(|v| v.len()).max().unwrap();

        Self {
            n_data,
	    n_total,
            m,
            n_edges: edge_to_vn.len(),
            cn_edges_data,
	    cn_edges_total,
            vn_edges,
            cn_max_deg,
            vn_max_deg,
            edge_to_vn,
	    _mode: PhantomData,
        }
    }
}

pub fn build_classical_graph(h_csr: &CsMat<u8>)
			     -> (Vec<Vec<usize>>, Vec<Vec<usize>>, Vec<Vec<usize>>, Vec<usize>) {

    let rows = h_csr.rows();
    let cols = h_csr.cols();

    let mut cn_edges_data = vec![Vec::new(); rows];
    let mut vn_edges = vec![Vec::new(); cols];
    let mut edge_to_vn = Vec::new();
    
    let binding = h_csr.indptr();
    let indptr  = binding.as_slice().unwrap();
    let indices = h_csr.indices();

    let mut edge_id = 0;

    for row in 0..rows {
        for idx in indptr[row]..indptr[row + 1] {
            let col = indices[idx];

            cn_edges_data[row].push(edge_id);
            vn_edges[col].push(edge_id);

	    edge_to_vn.push(col);
	    
            edge_id += 1;
        }
    }
    let cn_edges_total = cn_edges_data.clone(); // No graph extension

    (cn_edges_data, cn_edges_total, vn_edges, edge_to_vn)
}

pub fn build_quantum_graph(h_csr: &CsMat<u8>)
			   -> (Vec<Vec<usize>>, Vec<Vec<usize>>, Vec<Vec<usize>>, Vec<usize>) {

    let rows = h_csr.rows();
    let cols = h_csr.cols();

    let mut cn_edges_data  = vec![Vec::new(); rows];
    let mut cn_edges_total = vec![Vec::new(); rows];
    let mut vn_edges = vec![Vec::new(); cols + rows];
    let mut edge_to_vn = Vec::new();

    let binding = h_csr.indptr();
    let indptr  = binding.as_slice().unwrap();
    let indices = h_csr.indices();

    let mut edge_id = 0;

    // --- 1. Normal H edges (data qubits) ---
    for row in 0..rows {
        for idx in indptr[row]..indptr[row + 1] {
            let col = indices[idx];

            cn_edges_data[row].push(edge_id);
	    cn_edges_total[row].push(edge_id);
            vn_edges[col].push(edge_id);
            edge_to_vn.push(col);

            edge_id += 1;
        }
    }

    // --- 2. Identity edges (syndrome error variables) ---
    for j in 0..rows {
        let vn = cols + j; // e_j variable

        cn_edges_total[j].push(edge_id);
        vn_edges[vn].push(edge_id);
        edge_to_vn.push(vn);

        edge_id += 1;
    }

    (cn_edges_data, cn_edges_total, vn_edges, edge_to_vn)
}
