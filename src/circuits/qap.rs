use crate::math::{field_element::FieldElement as FE, polynomial::Polynomial};

#[derive(Clone, Debug, PartialEq, Eq)]
/// QAP Representation of the circuits
pub struct Qap {
    pub v: Vec<Polynomial>,
    pub w: Vec<Polynomial>,
    pub y: Vec<Polynomial>,
    pub target: Polynomial,
    pub number_of_inputs: usize,
    pub number_of_outputs: usize,
}

#[derive(Debug, PartialEq, Eq)]
pub enum CreationError {
    PolynomialVectorsSizeMismatch,
}

impl Qap {
    /// Creates a new QAP
    /// This expects vectors to be organized like:
    /// v0,w0,y0
    /// inputs associated v,w,y polynomials
    /// mid associated polynomials
    /// outputs associated v,w,y polynomials
    pub fn new(
        v: Vec<Polynomial>,
        w: Vec<Polynomial>,
        y: Vec<Polynomial>,
        target: Polynomial,
        number_of_inputs: usize,
        number_of_outputs: usize,
    ) -> Result<Self, CreationError> {
        // TO DO: Check if the amount of inputs and outputs matches the polynomials
        if v.len() != w.len() || v.len() != y.len() || w.len() != y.len() {
            Err(CreationError::PolynomialVectorsSizeMismatch)
        } else {
            Ok(Self {
                v,
                w,
                y,
                target,
                number_of_inputs,
                number_of_outputs,
            })
        }
    }

    pub fn v_mid(&'_ self) -> &[Polynomial] {
        &self.v[self.number_of_inputs + 1..(self.v.len() - self.number_of_outputs)]
    }

    pub fn w_mid(&'_ self) -> &[Polynomial] {
        &self.w[self.number_of_inputs + 1..(self.w.len() - self.number_of_outputs)]
    }

    pub fn y_mid(&'_ self) -> &[Polynomial] {
        &self.y[self.number_of_inputs + 1..(self.y.len() - self.number_of_outputs)]
    }
    
    pub fn v_input(&'_ self) -> &[Polynomial]{
        &self.y[1..self.number_of_inputs + 1]
    }

    pub fn w_input(&'_ self) -> &[Polynomial]{
        &self.y[1..self.number_of_inputs + 1]
    }

    pub fn y_input(&'_ self) -> &[Polynomial]{
        &self.y[1..self.number_of_inputs + 1]
    }

    pub fn v0(&'_ self) -> &Polynomial{
        &self.v[0]
    }

    pub fn w0(&'_ self) -> &Polynomial{
        &self.w[0]
    }

    pub fn y0(&'_ self) -> &Polynomial{
        &self.y[0]
    }

    pub fn v_output(&'_ self) -> &[Polynomial]{
        &self.y[(self.v.len() - self.number_of_outputs)..]
    }
    pub fn w_output(&'_ self) -> &[Polynomial]{
        &self.w[(self.w.len() - self.number_of_outputs)..]
    }

    pub fn y_output(&'_ self) -> &[Polynomial]{
        &self.y[(self.y.len() - self.number_of_outputs)..]
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

        Self::new(vs.to_vec(), ws.to_vec(), ys.to_vec(), t, 4, 1).unwrap()
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
            Qap::new(v, u, w, t, 2, 1)
        );
    }

    #[test]
    fn test_circuit_v_w_y_have_7_elements() {
        let test_circuit = Qap::new_test_circuit();
        assert_eq!(test_circuit.v.len(), 7);
        assert_eq!(test_circuit.w.len(), 7);
        assert_eq!(test_circuit.y.len(), 7);
    }

    //_mid polynomials of test circuit contains only one polynomial
    #[test]
    fn v_mid_test_circuit_on_r6_is_0() {
        let test_circuit = Qap::new_test_circuit();
        let r6 = FE::new(2).unwrap();
        assert_eq!(test_circuit.y_mid()[0].evaluate(r6), FE::zero());
    }

    #[test]
    fn w_mid_test_circuit_has_one_element() {
        let test_circuit = Qap::new_test_circuit();
        print!("mid w: {:?}", test_circuit.w_mid());
        assert_eq!(test_circuit.v_mid().len(), 1);
    }

    #[test]
    fn w_mid_test_circuit_on_r5_is_0() {
        let test_circuit = Qap::new_test_circuit();
        let r5 = FE::new(1).unwrap();
        assert_eq!(test_circuit.w_mid()[0].evaluate(r5), FE::zero());
    }

    #[test]
    fn w_mid_test_circuit_on_r6_is_1() {
        let test_circuit = Qap::new_test_circuit();
        let r6 = FE::new(2).unwrap();
        assert_eq!(test_circuit.w_mid()[0].evaluate(r6), FE::one());
    }

    #[test]
    fn y_mid_test_circuit_on_r5_is_1() {
        let test_circuit = Qap::new_test_circuit();
        let r5 = FE::new(1).unwrap();
        assert_eq!(test_circuit.y_mid()[0].evaluate(r5), FE::one());
    }

    #[test]
    fn y_mid_test_circuit_on_r6_is_0() {
        let test_circuit = Qap::new_test_circuit();
        let r6 = FE::new(2).unwrap();
        assert_eq!(test_circuit.y_mid()[0].evaluate(r6), FE::zero());
    }

    #[test]
    fn v_input_test_circuit_has_length_4() {
        let test_circuit = Qap::new_test_circuit();
        assert_eq!(test_circuit.v_input().len(), 4);
    }

    #[test]
    fn w_input_test_circuit_has_length_4() {
        let test_circuit = Qap::new_test_circuit();
        assert_eq!(test_circuit.w_input().len(), 4);
    }
    #[test]
    fn y_input_test_circuit_has_length_4() {
        let test_circuit = Qap::new_test_circuit();
        assert_eq!(test_circuit.y_input().len(), 4);
    }

    #[test]
    fn v_output_test_circuit_has_length_1() {
        let test_circuit = Qap::new_test_circuit();
        assert_eq!(test_circuit.v_output().len(), 1);
    }

    #[test]
    fn w_output_test_circuit_has_length_1() {
        let test_circuit = Qap::new_test_circuit();
        assert_eq!(test_circuit.w_output().len(), 1);
    }

    #[test]
    fn y_output_test_circuit_has_length_1() {
        let test_circuit = Qap::new_test_circuit();
        assert_eq!(test_circuit.y_output().len(), 1);
    }
}
