use crate::mode::{Quantum};
use crate::graph::Graph;
use crate::channel::{QuantumBscChannel};
use crate::random_ldpc::gen_ldpc;
use crate::decoder::Decoder;

pub fn run() {
    // Create regular random code, set simulation parameters
    let runs = 100;
    let n:usize  = 2000;
    let dv:usize = 3;
    let dc:usize = 6;
    let h_csr = gen_ldpc(n, dv, dc);

    // qLDPC joint decoding: test with QuantumBSC channel
    // p_error:    Pauli X-Error (bit flip)
    // p_syn_flip: Syndrome measurement error (syndrome bit flip)
    //
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

    for _ in 0..runs {
        match dec.apply_channel(&channel, &tx) {
            Ok(_) => (),
            Err(e) => println!("error: {}", e),
        }
        dec.decode();
        println!("{}", dec.result);
    }
}
