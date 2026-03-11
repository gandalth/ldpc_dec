mod decoder;
use decoder::Decoder;

mod random_ldpc;
use random_ldpc::gen_ldpc;

use rand::rng;
use rand_distr::{Distribution, StandardNormal};

fn main() {

    // Create regular random code
    let n:usize  = 2000;
    let dv:usize = 3;
    let dc:usize = 6;
    let h_csr = gen_ldpc(n, dv, dc);

    // Initialize decoder and provide basic info
    let mut dec = Decoder::new(h_csr, vec![]);
    dec.info();

    // Create sample received vector (AWGN output)
    let mut rng = rng();
    let sigma = 0.8;

    // recv: encoded, bpsk-mapped, AWGN-noise : 2 * mod(x*G, 2) - 1 + noise
    // sigma: std-dev of AWGN noise
    // For quick check, use all-zeros codeword
    let mut recv = vec![-1.0; n];
    for ri in &mut recv {
	let noise: f32 = StandardNormal.sample(&mut rng);
	*ri += sigma * noise;
    }


    match dec.decode(&recv, sigma) {
	Ok(_) => (),
	Err(e) => println!("error: {}", e),
    }
}




