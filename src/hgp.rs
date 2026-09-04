
use crate::matrix::{hcat, kronecker, print_matrix};
use crate::random_ldpc::{gen_ldpc};
use sprs::{CsMat};

// Create a random code and derive a hypergraph code from it

pub fn hgp_from_random(n_seed: usize, dv: usize, dc: usize)
                       -> (CsMat<u8>, CsMat<u8>) {
    let h        = gen_ldpc(n_seed, dv, dc);
    let (hx, hz) = hgp(&h, &h);

    (hx, hz)
}

// Calculate hypergraph product code from two matrices

pub fn hgp(h1: &CsMat<u8>, h2: &CsMat<u8>,) -> (CsMat<u8>, CsMat<u8>) {

    let (m1, n1) = h1.shape();
    let (m2, n2) = h2.shape();

    let i_n1 = CsMat::<u8>::eye(n1);
    let i_n2 = CsMat::<u8>::eye(n2);
    let i_m1 = CsMat::<u8>::eye(m1);
    let i_m2 = CsMat::<u8>::eye(m2);

    let h1_t = h1.transpose_view().to_owned().to_csr();
    let h2_t = h2.transpose_view().to_owned().to_csr();

    // Hx = [ H1 \kron I_n2 | I_m1 \kron H2^T ]

    let hx_left  = kronecker(h1, &i_n2);
    let hx_right = kronecker(&i_m1, &h2_t);

    let hx = hcat(&hx_left, &hx_right);

    // Hz = [ I_n1 \kron H2 | H1^T \kron I_m2 ]

    let hz_left  = kronecker(&i_n1, h2);
    let hz_right = kronecker(&h1_t, &i_m2);

    let hz = hcat(&hz_left, &hz_right);

    (hx, hz)
}
