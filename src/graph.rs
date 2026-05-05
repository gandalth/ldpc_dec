use sprs::CsMat;

pub struct Graph {
    pub n:          usize,
    pub m:          usize,
    pub n_edges:    usize,
    pub cn_edges:   Vec<Vec<usize>>,
    pub vn_edges:   Vec<Vec<usize>>,
    pub cn_max_deg: usize,
    pub vn_max_deg: usize,
    pub edge_to_vn: Vec<usize>
}


impl Graph {
    // Constructor
    pub fn new(h: CsMat<u8>) -> Self {

	let (m, n)  = h.shape();
	let n_edges = h.nnz();

	// Build the graph and the edges
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


pub fn build_graph(h_csr: &CsMat<u8>)
		   -> (Vec<Vec<usize>>, Vec<Vec<usize>>, Vec<usize>) {

    let rows = h_csr.rows();
    let cols = h_csr.cols();

    let mut cn_edges = vec![Vec::new(); rows];
    let mut vn_edges = vec![Vec::new(); cols];
    let mut edge_to_vn = Vec::new();
    
    let binding = h_csr.indptr();
    let indptr = binding.as_slice().unwrap();
    let indices = h_csr.indices();

    let mut edge_id = 0;

    for row in 0..rows {
        for idx in indptr[row]..indptr[row + 1] {
            let col = indices[idx];

            cn_edges[row].push(edge_id);
            vn_edges[col].push(edge_id);

	    edge_to_vn.push(col);
	    
            edge_id += 1;
        }
    }
    (cn_edges, vn_edges, edge_to_vn)
}
