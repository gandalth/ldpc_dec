mod decoder;
use decoder::{Decoder, OpMode};

mod graph;
use graph::Graph;

mod random_ldpc;
use random_ldpc::gen_ldpc;

mod channel;
use channel::{AwgnChannel, QuantumBscChannel};

mod node_math;

fn main() {

    // Create regular random code
    let n:usize  = 2000;
    let dv:usize = 3;
    let dc:usize = 6;
    let h_csr = gen_ldpc(n, dv, dc);
    let mode = OpMode::Classic;

    // Create graph
    let graph = Graph::new(h_csr);

    // Initialize decoder and provide basic info
    let mut dec = Decoder::new(&graph, mode, vec![]);
    dec.info();

    // tx: all-zeros, encoded, bpsk-mapped: 2 * mod(x*G, 2) - 1
    // sigma: std-dev of AWGN noise
    let sigma = 0.8;
    let runs = 100;
    let tx = vec![0u8; n];

    // let channel = AwgnChannel { sigma };
    let channel = QuantumBscChannel {
	graph: &graph,
	p_error: 0.005,
	p_syn_flip: 0.001,
    };

    for _ in 0..runs {

	match dec.apply_channel(&channel, &tx) {
	    Ok(_) => (),
	    Err(e) => println!("error: {}", e),
	}

	dec.decode();
	println!("{}", dec.result);
    }
}




