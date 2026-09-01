use crate::mode::{Quantum};
use crate::graph::Graph;
use crate::channel::{QuantumBscChannel};
use crate::css::{CssCode, CssGraphs, css_3qubit, css_steane_code};
use crate::decoder::Decoder;

pub fn run() {
    // CSS decoding: test with CSS-BSC channel
    // p_x_error:    Pauli X-Error (bit flip)
    // p_x_syn_flip: Syndrome measurement error (on measuring x errors)
    // p_z_error:    Pauli Z-Error (phase flip)
    // p_z_syn_flip: Syndrome measurement error (on measuring z errors)

    let runs = 100;
    let n:usize  = 7;
    let (hx, hz) = css_steane_code();

    let css = CssCode::new(hx, hz)
        .expect("CSS validity check failed.");
    let graphs = CssGraphs::new(&css);

    let tx = vec![0; n];

    // Create two independent QuantumBSC Channels

    let x_channel = QuantumBscChannel {
        graph: &graphs.z, // Hz detects X errors
        p_error: 0.1,
        p_syn_flip: 0.0,
    };
        
    let z_channel = QuantumBscChannel {
        graph: &graphs.x, // Hx detects Z errors
        p_error: 0.0,
        p_syn_flip: 0.0,
    };

    let mut dec_z = Decoder::<Quantum>::new(&graphs.x, vec![]);
    let mut dec_x = Decoder::<Quantum>::new(&graphs.z, vec![]);

    dec_z.info();
    dec_x.info();

    for _ in 0..runs {
        match dec_x.apply_channel(&x_channel, &tx)  {
            Ok(_) => (),
            Err(e) => println!("X-Decoder Error: {}", e),
        }
        match dec_z.apply_channel(&z_channel, &tx)  {
            Ok(_) => (),
            Err(e) => println!("Z-Decoder Error: {}", e),
        }

        dec_x.decode();
        println!("X-Decoder Result: {}", dec_x.result);
        dec_z.decode();
        println!("Z-Decoder Result: {}", dec_z.result);
    } 
}
