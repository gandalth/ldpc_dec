use crate::mode::{Classic};
use crate::graph::Graph;
use crate::channel::{AwgnChannel};
use crate::random_ldpc::gen_ldpc;
use crate::decoder::Decoder;

pub fn run() {
    // Create regular random code, set simulation parameters
    let runs = 100;
    let n:usize  = 2000;
    let dv:usize = 3;
    let dc:usize = 6;
    let h_csr = gen_ldpc(n, dv, dc);

    // Classic LDPC decoder: test with simple AWGN channel
    // sigma: std-dev of AWGN noise
    
    // Create graph
    let graph = Graph::<Classic>::new(h_csr);

    // tx: all-zeros, encoded, bpsk-mapped: 2 * mod(x*G, 2) - 1
    let tx = vec![-1.0; n];
    let channel = AwgnChannel { sigma : 0.8 };

    // Initialize decoder
    let mut dec = Decoder::<Classic>::new(&graph, vec![]);
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
