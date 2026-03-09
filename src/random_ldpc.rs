use sprs::{CsMat, TriMat};
use rand::seq::SliceRandom;
use rand::rng;

pub fn gen_ldpc(n: usize, dv: usize, dc: usize) -> CsMat<u8> {

    assert!(n * dv % dc == 0);

    let m = n * dv / dc;
    let edges = n * dv;

    let mut tri = TriMat::<u8>::with_capacity((m, n), edges);

    // create socket lists
    let mut vn_sockets = Vec::with_capacity(edges);
    let mut cn_sockets = Vec::with_capacity(edges);

    for v in 0..n {
        for _ in 0..dv {
            vn_sockets.push(v);
        }
    }

    for c in 0..m {
        for _ in 0..dc {
            cn_sockets.push(c);
        }
    }

    // random permutation (Gallager style)
    let mut rng = rng();
    cn_sockets.shuffle(&mut rng);

    for i in 0..edges {
        let v = vn_sockets[i];
        let c = cn_sockets[i];
        tri.add_triplet(c, v, 1);
    }

    tri.to_csr()
}
