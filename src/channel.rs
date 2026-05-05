use rand_distr::{Distribution, StandardNormal};

use crate::decoder::DecoderState;
use crate::node_math::hard_decision;

pub trait Channel {
    fn apply(&self, tx: &[f32], state: &mut DecoderState) -> Result<(), String>;
}

pub struct AwgnChannel {
    pub sigma: f32,
}

impl Channel for AwgnChannel {
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
        state.soft_syndrome.fill(1.0);
	state.hard_syndrome = hard_decision(&state.soft_syndrome);
	Ok(())
    }
}
