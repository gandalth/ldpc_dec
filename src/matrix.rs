use sprs::{CsMat, TriMat};

pub fn kronecker(a: &CsMat<u8>, b: &CsMat<u8>) -> CsMat<u8> {
    // Note: both matrices have to be in CSR representation.
    // Uncomment below if that is not the case.
    // let a = a.to_csr();
    // let b = b.to_csr();

    let (arows, acols) = a.shape();
    let (brows, bcols) = b.shape();

    let mut result =
        TriMat::<u8>::new((arows * brows, acols * bcols));
    
    for (ai, arow) in a.outer_iterator().enumerate() {
        for (aj, &aval) in arow.iter() {
            if aval == 0 {
                continue;
            }

            for (bi, brow) in b.outer_iterator().enumerate() {
                for (bj, &bval) in brow.iter() {
                    if bval == 0 {
                        continue;
                    }

                    let row = ai * brows + bi;
                    let col = aj * bcols + bj;

                    assert!(
                        row < arows * brows,
                        "row={} but result has {} rows",
                        row,
                        arows * brows
                    );

                    assert!(
                        col < acols * bcols,
                        "col={} but result has {} cols",
                        col,
                        acols * bcols
                    );
                    
                    result.add_triplet(row, col, aval * bval);
                }
            }
        }
    }

    result.to_csr()
}

pub fn hcat(a: &CsMat<u8>, b: &CsMat<u8>,) -> CsMat<u8> {

    assert_eq!(a.rows(), b.rows(),
        "hcat: number of rows not equal: {}, {}", a.rows(), b.rows()
    );

    // Note: both matrices have to be in CSR representation.
    // Uncomment below if that is not the case.
    // let a = a.to_csr();
    // let b = b.to_csr();
        
    let rows = a.rows();
    let cols = a.cols() + b.cols();

    let mut result = TriMat::<u8>::new((rows, cols));

    for (i, row) in a.outer_iterator().enumerate() {
        for (j, &value) in row.iter() {
            result.add_triplet(i, j, value);
        }
    }

    let offset = a.cols();

    for (i, row) in b.outer_iterator().enumerate() {
        for (j, &value) in row.iter() {
            result.add_triplet(i, offset + j, value);
        }
    }

    result.to_csr()
}

pub fn print_matrix(name: &str, m: &CsMat<u8>) {
    println!("{name} ({} × {}):", m.rows(), m.cols());

    for i in 0..m.rows() {
        for j in 0..m.cols() {
            print!("{} ", m.get(i, j).copied().unwrap_or(0));
        }
        println!();
    }
    println!();
}
