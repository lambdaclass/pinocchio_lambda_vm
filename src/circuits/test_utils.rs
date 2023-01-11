use super::qap::Qap;
use crate::math::{field_element::FieldElement as FE, polynomial::Polynomial};

// r5 and r6 are exposed to help testing
pub fn test_qap_r5() -> FE {
    FE::new(0)
}

pub fn test_qap_r6() -> FE {
    FE::new(1)
}

/// This is a solver for the test qap
/// Inputs: c1,c2,c3,c4 circuit inputs
/// Outputs: c5 intermediate result, c6 result
pub fn test_qap_solver(inputs: [FE; 4]) -> (FE, FE) {
    let c5 = inputs[2] * inputs[3];
    let c6 = (inputs[0] + inputs[1]) * c5;
    (c5, c6)
}

/// Test qap based on pinocchios paper example
pub fn new_test_qap() -> Qap {
    let r5: FE = test_qap_r5();
    let r6: FE = test_qap_r6();

    let t: Polynomial =
        Polynomial::new(vec![-r5, FE::new(1)]) * Polynomial::new(vec![-r6, FE::new(1)]);

    let vs = &[
        // v0 is 0 for everypoint for the circuit, since it has no constants
        // in the paper they don't write it
        Polynomial::interpolate(&[r5, r6], &[FE::new(0), FE::new(0)]),
        // v1..v6 are the ones explicitly written in the paper
        Polynomial::interpolate(&[r5, r6], &[FE::new(0), FE::new(1)]),
        Polynomial::interpolate(&[r5, r6], &[FE::new(0), FE::new(1)]),
        Polynomial::interpolate(&[r5, r6], &[FE::new(1), FE::new(0)]),
        Polynomial::interpolate(&[r5, r6], &[FE::new(0), FE::new(0)]),
        Polynomial::interpolate(&[r5, r6], &[FE::new(0), FE::new(0)]),
        Polynomial::interpolate(&[r5, r6], &[FE::new(0), FE::new(0)]),
    ];

    let ws = &[
        //w0
        Polynomial::interpolate(&[r5, r6], &[FE::new(0), FE::new(0)]),
        //w1
        Polynomial::interpolate(&[r5, r6], &[FE::new(0), FE::new(0)]),
        Polynomial::interpolate(&[r5, r6], &[FE::new(0), FE::new(0)]),
        Polynomial::interpolate(&[r5, r6], &[FE::new(0), FE::new(0)]),
        Polynomial::interpolate(&[r5, r6], &[FE::new(1), FE::new(0)]),
        Polynomial::interpolate(&[r5, r6], &[FE::new(0), FE::new(1)]),
        Polynomial::interpolate(&[r5, r6], &[FE::new(0), FE::new(0)]),
    ];

    let ys = &[
        //y0
        Polynomial::interpolate(&[r5, r6], &[FE::new(0), FE::new(0)]),
        //y1
        Polynomial::interpolate(&[r5, r6], &[FE::new(0), FE::new(0)]),
        Polynomial::interpolate(&[r5, r6], &[FE::new(0), FE::new(0)]),
        Polynomial::interpolate(&[r5, r6], &[FE::new(0), FE::new(0)]),
        Polynomial::interpolate(&[r5, r6], &[FE::new(0), FE::new(0)]),
        Polynomial::interpolate(&[r5, r6], &[FE::new(1), FE::new(0)]),
        Polynomial::interpolate(&[r5, r6], &[FE::new(0), FE::new(1)]),
    ];

    Qap::new(vs.to_vec(), ws.to_vec(), ys.to_vec(), t, 4, 1).unwrap()
}
