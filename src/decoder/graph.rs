use sprs::CsMat;

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
