pub fn normalized_mult_exc_one(f0: &[f32],
			       prefix_f0: &mut [f32], suffix_f0: &mut [f32],
			       prefix_f1: &mut [f32], suffix_f1: &mut [f32],
			       result: &mut [f32]) {

    let n = f0.len();

    prefix_f0.fill(1.0);
    suffix_f0.fill(1.0);
    prefix_f1.fill(1.0);
    suffix_f1.fill(1.0);

    for i in 1..n {
        prefix_f0[i] = prefix_f0[i - 1] * f0[i - 1];
	prefix_f1[i] = prefix_f1[i - 1] * (1.0 - f0[i - 1]);
    }

    for i in (0..n - 1).rev() {
        suffix_f0[i] = suffix_f0[i + 1] * f0[i + 1];
	suffix_f1[i] = suffix_f1[i + 1] * (1.0 - f0[i + 1]);
    }

    for i in 0..n {
	let num = prefix_f0[i] * suffix_f0[i];
	let den = num + prefix_f1[i] * suffix_f1[i];
	if den == 0.0 {
		result[i] = 0.5
	    } else {
		result[i] = num / den
	}
    }
}

pub fn gallager_prod_exc_one(f0: &[f32],
			     prefix_f0: &mut [f32], suffix_f0: &mut [f32],
			     result: &mut [f32]) {
    // Calculate Gallager product over input slice leaving out exactly one
    // input element in each variation (to receive extrinsic information)

    let n = f0.len();

    prefix_f0.fill(1.0);
    suffix_f0.fill(1.0);

    for i in 1..n {
	let x_prev = 2.0 * f0[i - 1] - 1.0;
	prefix_f0[i]  = prefix_f0[i - 1] * x_prev;
    }

    for i in (0..n - 1).rev() {
	let x_next = 2.0 * f0[i + 1] - 1.0;
        suffix_f0[i]  = suffix_f0[i + 1] * x_next;
    }

    for i in 0..n {
	result[i] = 0.5 * (prefix_f0[i] * suffix_f0[i]) + 0.5;
    }
}

pub fn normalized_mult(f0: &[f32]) -> f32 {
   
    let f1: Vec<f32> = f0
	.iter()
	.map(|&r| 1.0 - r)
	.collect();

    let prod_f0: f32 = f0.iter().product();
    let prod_f1: f32 = f1.iter().product();
    let num = prod_f0;
    let den = prod_f0 + prod_f1;
    
    if den == 0.0 {
	0.5
    } else {
	num / den
    }
 
}

pub fn hard_decision(p0: &[f32]) -> Vec<u8> {
    p0.iter().map(|&p| (p < 0.5) as u8).collect()
}

fn _gallager_prod(f0: &[f32]) -> f32 {
    let p0: f32 = f0
	.iter()
	.map(|&f0| 2.0 * f0 - 1.0)
	.product();
    return 0.5 * p0 + 0.5	     
}

