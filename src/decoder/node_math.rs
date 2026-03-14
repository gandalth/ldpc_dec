pub fn normalized_mult_exc_one(f0: &[f32]) -> Vec<f32> {

    let n = f0.len();
    if n == 0 {
	return Vec::new();
    }
    
    let mut prefix_f0 = vec![1.0; n];
    let mut suffix_f0 = vec![1.0; n];
    let mut prefix_f1 = vec![1.0; n];
    let mut suffix_f1 = vec![1.0; n];
				    
    for i in 1..n {
        prefix_f0[i] = prefix_f0[i - 1] * f0[i - 1];
	prefix_f1[i] = prefix_f1[i - 1] * (1.0 - f0[i - 1]);
    }
    
    for i in (0..n - 1).rev() {
        suffix_f0[i] = suffix_f0[i + 1] * f0[i + 1];
	suffix_f1[i] = suffix_f1[i + 1] * (1.0 - f0[i + 1]);
    }
    
    (0..n)
        .map(|i| {
	    let num = prefix_f0[i] * suffix_f0[i];
	    let den = num + prefix_f1[i] * suffix_f1[i];
	    if den == 0.0 {
		0.5
	    } else {
		num / den
	    }
	})
        .collect()
}

pub fn gallager_prod_exc_one(f0: &[f32]) -> Vec<f32> {
    // Calculate Gallager product over input slice leaving out exactly one
    // input element in each variation (to receive extrinsic information)

    let n = f0.len();
    if n == 0 {
	return Vec::new();
    }
    
    let mut prefix = vec![1.0; n];
    let mut suffix = vec![1.0; n];
    let mut result = vec![0.0; n];

    for i in 1..n {
	let x_prev = 2.0 * f0[i - 1] - 1.0;
	prefix[i]  = prefix[i - 1] * x_prev;
    }

    for i in (0..n - 1).rev() {
	let x_next = 2.0 * f0[i + 1] - 1.0;
        suffix[i]  = suffix[i + 1] * x_next;
    }

    for i in 0..n {
	result[i] = 0.5 * (prefix[i] * suffix[i]) + 0.5;
    }
    result
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

