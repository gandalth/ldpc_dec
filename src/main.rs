mod decoder;
use decoder::{Decoder, OpMode};

mod random_ldpc;
use random_ldpc::gen_ldpc;

mod channel;
use channel::{AwgnChannel};

mod node_math;

fn main() {

    // Create regular random code
    let n:usize  = 2000;
    let dv:usize = 3;
    let dc:usize = 6;
    let h_csr = gen_ldpc(n, dv, dc);
    let mode = OpMode::Classic;

    // Initialize decoder and provide basic info
    let mut dec = Decoder::new(h_csr, mode, vec![]);
    dec.info();

    // Create sample received vector (AWGN output)
    let sigma = 0.8;
    let runs = 100;

    // tx: all-zeros, encoded, bpsk-mapped: 2 * mod(x*G, 2) - 1
    // sigma: std-dev of AWGN noise
    let tx = vec![-1.0; n];

    let channel = AwgnChannel { sigma };

    for _ in 0..runs {

	match dec.apply_channel(&channel, &tx) {
	    Ok(_) => (),
	    Err(e) => println!("error: {}", e),
	}

	dec.decode();
	println!("{}", dec.result);
    }
}




