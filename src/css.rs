use sprs::{CsMat,TriMat};
use crate::graph::Graph;
use crate::mode::Quantum;

pub struct CssCode {
    pub hx: CsMat<u8>,
    pub hz: CsMat<u8>,
}

impl CssCode {
    // Constructor: check Hx * Hz^T = 0 (mod 2) and build CSS.
    pub fn new(hx: CsMat<u8>, hz: CsMat<u8>) -> Result<Self, String> {

        check_css(&hx, &hz)?;

        Ok(Self { hx, hz })
    }
}

pub struct CssGraphs {
    pub x: Graph<Quantum>,
    pub z: Graph<Quantum>,
}

impl CssGraphs {
    pub fn new(code: &CssCode) -> Self {

        let x = Graph::<Quantum>::new(code.hx.clone());
        let z = Graph::<Quantum>::new(code.hz.clone());

        Self {
            x,
            z,
        }
    }
}


// 3-qubit CSS code
pub fn css_3qubit() -> (CsMat<u8>, CsMat<u8>) {

    let mut hx = TriMat::<u8>::new((2, 3));

    // row 0: 1 1 0
    hx.add_triplet(0, 0, 1);
    hx.add_triplet(0, 1, 1);

    // row 1: 0 1 1
    hx.add_triplet(1, 1, 1);
    hx.add_triplet(1, 2, 1);

    let hx = hx.to_csr();

    let mut hz = TriMat::<u8>::new((1, 3));

    // row 0: 1 1 1
    hz.add_triplet(0, 0, 1);
    hz.add_triplet(0, 1, 1);
    hz.add_triplet(0, 2, 1);

    let hz = hz.to_csr();

    // Return pair of parity-check matrices
    (hx, hz)
}

pub fn css_steane_code() -> (CsMat<u8>, CsMat<u8>) {

    let mut h = TriMat::<u8>::new((3, 7));

    // Row 0: 1 0 0 1 0 1 1
    h.add_triplet(0, 0, 1);
    h.add_triplet(0, 3, 1);
    h.add_triplet(0, 5, 1);
    h.add_triplet(0, 6, 1);

    // Row 1: 0 1 0 1 1 0 1
    h.add_triplet(1, 1, 1);
    h.add_triplet(1, 3, 1);
    h.add_triplet(1, 4, 1);
    h.add_triplet(1, 6, 1);

    // Row 2: 0 0 1 0 1 1 1
    h.add_triplet(2, 2, 1);
    h.add_triplet(2, 4, 1);
    h.add_triplet(2, 5, 1);
    h.add_triplet(2, 6, 1);

    let h = h.to_csr();

    // Return pair of parity-check matrices
    (h.clone(), h)
}

// check_css: check orthogonality constraint for CSS codes
pub fn check_css(hx: &CsMat<u8>, hz: &CsMat<u8>) -> Result<(), String> {
    if hx.cols() != hz.cols() {
        return Err(format!(
            "Hx and Hz must have the same number of columns (got {} and {})",
            hx.cols(),
            hz.cols()
        ));
    }

    for (ix, row_x) in hx.outer_iterator().enumerate() {

        for (iz, row_z) in hz.outer_iterator().enumerate() {

            let mut parity = 0u8;

            let mut px = 0;
            let mut pz = 0;

            let idx_x = row_x.indices();
            let idx_z = row_z.indices();

            while px < idx_x.len() && pz < idx_z.len() {

                if idx_x[px] == idx_z[pz] {
                    parity ^= 1;
                    px += 1;
                    pz += 1;
                }
                else if idx_x[px] < idx_z[pz] {
                    px += 1;
                }
                else {
                    pz += 1;
                }
            }

            if parity != 0 {
                return Err(format!(
                    "CSS condition violated:\
                     row {} of Hx and row {} of Hz anticommute", ix, iz));
            }
        }
    }

    Ok(())
}
