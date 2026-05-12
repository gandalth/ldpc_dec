use rand_distr::{Distribution, StandardNormal};

use crate::decoder::DecoderState;
use crate::graph::Graph;

use crate::mode::{Mode, Classic, Quantum};

use rand::Rng;

pub trait Channel<M: Mode> {
    type Tx;
    fn apply(&self, tx: &[Self::Tx], state: &mut DecoderState)
	     -> Result<(), String>;
}

pub struct AwgnChannel {
    pub sigma: f32,
}

impl Channel<Classic> for AwgnChannel {
    type Tx = f32;
    fn apply(&self, tx: &[f32], state: &mut DecoderState)
	     -> Result<(), String> {

	if tx.len() !=  state.p0_aprio.len() {
            return Err(format!(
		"AWGN apply(): tx / decoder block length mismatch {} <-> {}",
		tx.len(), state.p0_aprio.len()));
	}
	let mut rng = rand::rng();
	let alpha = 2.0 / (self.sigma * self.sigma);


	for (i, &x) in tx.iter().enumerate() {
            let noise: f32 = StandardNormal.sample(&mut rng);
            let y = x + self.sigma * noise;
            state.p0_aprio[i] = 1.0 / (1.0 + (alpha * y).exp());
        }

        // all parity checks are even
        state.syndrome.fill(0);
	Ok(())
    }
}

pub struct QuantumBscChannel<'a> {
    pub graph: &'a Graph<Quantum>,
    pub p_error: f32,       // physical qubit flip rate
    pub p_syn_flip: f32,    // syndrome measurement noise
}

impl<'a> Channel<Quantum> for QuantumBscChannel<'a> {
    type Tx = u8;
    fn apply(&self, tx: &[u8], state: &mut DecoderState)
	     -> Result<(), String> {


	// The quantum BSC channel flips variables, calculates the correct
	// syndrome and applies bit-flips to the syndrome as well.

	let mut rng = rand::rng();

        let n = self.graph.n_data;
        let m = self.graph.m;

	if tx.len() != n {
            return Err("tx length mismatch".to_string());
        }

        let mut error = vec![0u8; n];

        for i in 0..n {
            if rng.random::<f32>() < self.p_error {
                error[i] = 1;
            }
        }

        let mut rx = vec![0u8; n];

        for i in 0..n {
            rx[i] = tx[i] ^ error[i];
        }

        let mut syndrome = vec![0u8; m];

        for (j, edges) in self.graph.cn_edges_data.iter().enumerate() {
            let mut parity = 0u8;

            for &e in edges {
                let v = self.graph.edge_to_vn[e];
                parity ^= error[v];
	    }

            syndrome[j] = parity;
        }

        for j in 0..m {
            if rng.random::<f32>() < self.p_syn_flip {
                syndrome[j] ^= 1;
            }
        }

	// Set priors for variable nodes
	for i in 0..n {
	    state.p0_aprio[i] = 1.0 - self.p_error;
	}

	// Store yndrome errors in additional variable nodes
	for j in 0..m {
	    state.p0_aprio[n + j] =
		if syndrome[j] == 0
	    {1.0 - self.p_syn_flip} else {self.p_syn_flip};
	}

	// all parity checks are even
        state.syndrome.fill(0); // Noisy syndrome in extension of p0_aprio 

	Ok(())
    }
}
