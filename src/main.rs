mod decoder;
use decoder::Decoder;

mod graph;
use graph::Graph;

mod random_ldpc;
use random_ldpc::gen_ldpc;

mod channel;
use channel::{AwgnChannel, QuantumBscChannel};

mod mode;
use mode::{Classic, Quantum};

mod node_math;

fn main() {

    
    // Create regular random code, set simulation parameters
    let n:usize  = 2000;
    let dv:usize = 3;
    let dc:usize = 6;
    let h_csr = gen_ldpc(n, dv, dc);
    let runs = 100;

    // ---
    // qLDPC joint decoding: test with QuantumBSC channel
    // p_error:    Pauli X-Error (bit flip)
    // p_syn_flip: Syndrome measurement error (syndrome bit flip)

    // Create graph
    let graph = Graph::<Quantum>::new(h_csr);

    let tx = vec![0; n];
    let channel = QuantumBscChannel {
	graph: &graph,
	p_error: 0.025,
	p_syn_flip: 0.005,
    };

    // Initialize decoder
    let mut dec = Decoder::<Quantum>::new(&graph, vec![]);
    // ---

    //// ---
    //// Classic LDPC decoder: test with simple AWGN channel
    //// sigma: std-dev of AWGN noise
    //
    //// Create graph
    //let graph = Graph::<Classic>::new(h_csr);
    //
    //// tx: all-zeros, encoded, bpsk-mapped: 2 * mod(x*G, 2) - 1
    // let tx = vec![-1.0; n];
    // let channel = AwgnChannel { sigma : 0.8 };
    //
    //// Initialize decoder
    //let mut dec = Decoder::<Classic>::new(&graph, vec![]);
    //// ---


    // Provide info of decoder and code
    dec.info();
    
    for _ in 0..runs {

	match dec.apply_channel(&channel, &tx) {
	    Ok(_) => (),
	    Err(e) => println!("error: {}", e),
	}

	dec.decode();
	println!("{}", dec.result);
    }
}




