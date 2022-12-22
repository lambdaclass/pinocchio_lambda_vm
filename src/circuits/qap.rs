use crate::math::{field_element::FieldElement as FE, polynomial::Polynomial};
#[derive(Debug, PartialEq, Eq)]
pub struct Qap {
    pub v: Vec<Polynomial>,
    pub w: Vec<Polynomial>,
    pub y: Vec<Polynomial>,
    pub target: Polynomial,
}

#[derive(Debug, PartialEq, Eq)]
pub enum CreationError {
    PolynomialVectorsSizeMismatch,
}

/// QAP Representation of the circuits

impl Qap {
    pub fn new(
        v: Vec<Polynomial>,
        w: Vec<Polynomial>,
        y: Vec<Polynomial>,
        target: Polynomial,
    ) -> Result<Self, CreationError> {
        if v.len() != w.len() || v.len() != y.len() || w.len() != y.len() {
            Err(CreationError::PolynomialVectorsSizeMismatch)
        } else {
            Ok(Self { v, w, y, target })
        }
    }

    pub fn new_test_circuit() -> Self {
        let r5: FE = FE::new(1).unwrap();
        let r6: FE = FE::new(2).unwrap();

        let t: Polynomial = Polynomial::new(vec![-r5, FE::new(1).unwrap()])
            * Polynomial::new(vec![-r6, FE::new(1).unwrap()]);

        let vs = &[
            // v0 is 0 for everypoint for the circuit, since it has no constants
            // in the paper they don't write it
            Polynomial::interpolate(&[r5, r6], &[FE::zero(), FE::zero()]),
            // v1..v6 are the ones explicitly written in the paper
            Polynomial::interpolate(&[r5, r6], &[FE::zero(), FE::one()]),
            Polynomial::interpolate(&[r5, r6], &[FE::zero(), FE::one()]),
            Polynomial::interpolate(&[r5, r6], &[FE::one(), FE::zero()]),
            Polynomial::interpolate(&[r5, r6], &[FE::zero(), FE::zero()]),
            Polynomial::interpolate(&[r5, r6], &[FE::zero(), FE::zero()]),
            Polynomial::interpolate(&[r5, r6], &[FE::zero(), FE::zero()]),
        ];

        let ws = &[
            //w0
            Polynomial::interpolate(&[r5, r6], &[FE::zero(), FE::zero()]),
            //w1
            Polynomial::interpolate(&[r5, r6], &[FE::zero(), FE::zero()]),
            Polynomial::interpolate(&[r5, r6], &[FE::zero(), FE::zero()]),
            Polynomial::interpolate(&[r5, r6], &[FE::zero(), FE::zero()]),
            Polynomial::interpolate(&[r5, r6], &[FE::one(), FE::zero()]),
            Polynomial::interpolate(&[r5, r6], &[FE::zero(), FE::one()]),
            Polynomial::interpolate(&[r5, r6], &[FE::zero(), FE::zero()]),
        ];

        let ys = &[
            //y0
            Polynomial::interpolate(&[r5, r6], &[FE::zero(), FE::zero()]),
            //y1
            Polynomial::interpolate(&[r5, r6], &[FE::zero(), FE::zero()]),
            Polynomial::interpolate(&[r5, r6], &[FE::zero(), FE::zero()]),
            Polynomial::interpolate(&[r5, r6], &[FE::zero(), FE::zero()]),
            Polynomial::interpolate(&[r5, r6], &[FE::zero(), FE::zero()]),
            Polynomial::interpolate(&[r5, r6], &[FE::one(), FE::zero()]),
            Polynomial::interpolate(&[r5, r6], &[FE::zero(), FE::one()]),
        ];

        Self::new(vs.to_vec(), ws.to_vec(), ys.to_vec(), t).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn qap_with_different_amount_of_polynomials_should_error() {
        let v = vec![
            Polynomial::new(vec![FE::new(1).unwrap(), FE::new(2).unwrap()]),
            Polynomial::new(vec![FE::new(2).unwrap(), FE::new(3).unwrap()]),
        ];

        let u = v.clone();
        let w = vec![Polynomial::new(vec![
            FE::new(1).unwrap(),
            FE::new(2).unwrap(),
        ])];
        let t = Polynomial::new(vec![FE::new(3).unwrap()]);
        assert_eq!(
            Err(CreationError::PolynomialVectorsSizeMismatch),
            Qap::new(v, u, w, t)
        );
    }

    #[test]
    fn test_circuit_v0_w0_y0_have_7_elements() {
        let test_circuit = Qap::new_test_circuit();
        assert_eq!(test_circuit.v.len(), 7);
        assert_eq!(test_circuit.w.len(), 7);
        assert_eq!(test_circuit.y.len(), 7);
    }
}
