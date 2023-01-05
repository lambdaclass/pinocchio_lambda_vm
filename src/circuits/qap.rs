use crate::math::{field_element::FieldElement, polynomial::Polynomial};

const ORDER: u128 = 13;
type FE = FieldElement<ORDER>;

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

    pub fn h_polynomial(&self, c: &[FE]) -> Polynomial {
        self.p_polynomial(c).div_with_ref(&self.target)
    }
    /// Receives C elements of a solution of the circuit
    /// Returns p polynomial
    // This along the polynomial execution should be migrated with a better
    // representation of the circuit
    pub fn p_polynomial(&self, c: &[FE]) -> Polynomial {
        let v: Polynomial = self.v[0].clone()
            + self.v[1..]
                .iter()
                .zip(c)
                .map(|(v, c)| v.mul_with_ref(&Polynomial::new_monomial(*c, 0)))
                .reduce(|x, y| x + y)
                .unwrap();

        let w: Polynomial = self.w[0].clone()
            + self.w[1..]
                .iter()
                .zip(c)
                .map(|(w, c)| w.mul_with_ref(&Polynomial::new_monomial(*c, 0)))
                .reduce(|x, y| x + y)
                .unwrap();

        let y: Polynomial = self.y[0].clone()
            + self.y[1..]
                .iter()
                .zip(c)
                .map(|(y, c)| y.mul_with_ref(&Polynomial::new_monomial(*c, 0)))
                .reduce(|x, y| x + y)
                .unwrap();

        v * w - y
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

    pub fn v_input(&'_ self) -> &[Polynomial] {
        &self.v[1..self.number_of_inputs + 1]
    }

    pub fn w_input(&'_ self) -> &[Polynomial] {
        &self.w[1..self.number_of_inputs + 1]
    }

    pub fn y_input(&'_ self) -> &[Polynomial] {
        &self.y[1..self.number_of_inputs + 1]
    }

    pub fn v0(&'_ self) -> &Polynomial {
        &self.v[0]
    }

    pub fn w0(&'_ self) -> &Polynomial {
        &self.w[0]
    }

    pub fn y0(&'_ self) -> &Polynomial {
        &self.y[0]
    }

    pub fn v_output(&'_ self) -> &[Polynomial] {
        &self.v[(self.v.len() - self.number_of_outputs)..]
    }
    pub fn w_output(&'_ self) -> &[Polynomial] {
        &self.w[(self.w.len() - self.number_of_outputs)..]
    }

    pub fn y_output(&'_ self) -> &[Polynomial] {
        &self.y[(self.y.len() - self.number_of_outputs)..]
    }

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
}

/// Test qap based on pinocchios paper example
pub fn new_test_qap() -> Qap {
    let r5: FE = Qap::test_qap_r5();
    let r6: FE = Qap::test_qap_r6();

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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn qap_with_different_amount_of_polynomials_should_error() {
        let v = vec![
            Polynomial::new(vec![FE::new(1), FE::new(2)]),
            Polynomial::new(vec![FE::new(2), FE::new(3)]),
        ];

        let u = v.clone();
        let w = vec![Polynomial::new(vec![FE::new(1), FE::new(2)])];
        let t = Polynomial::new(vec![FE::new(3)]);
        assert_eq!(
            Err(CreationError::PolynomialVectorsSizeMismatch),
            Qap::new(v, u, w, t, 2, 1)
        );
    }

    #[test]
    fn test_circuit_v_w_y_have_7_elements() {
        let test_circuit = new_test_qap();
        assert_eq!(test_circuit.v.len(), 7);
        assert_eq!(test_circuit.w.len(), 7);
        assert_eq!(test_circuit.y.len(), 7);
    }

    //_mid polynomials of test circuit contains only one polynomial
    #[test]
    fn v_mid_test_circuit_on_r6_is_0() {
        let test_circuit = new_test_qap();
        let r6 = Qap::test_qap_r6();
        assert_eq!(test_circuit.y_mid()[0].evaluate(r6), FE::new(0));
    }

    #[test]
    fn w_mid_test_circuit_has_one_element() {
        let test_circuit = new_test_qap();
        assert_eq!(test_circuit.v_mid().len(), 1);
    }

    #[test]
    fn w_mid_test_circuit_on_r5_is_0() {
        let test_circuit = new_test_qap();
        let r5 = Qap::test_qap_r5();
        assert_eq!(test_circuit.w_mid()[0].evaluate(r5), FE::new(0));
    }

    #[test]
    fn w_mid_test_circuit_on_r6_is_1() {
        let test_circuit = new_test_qap();
        let r6 = Qap::test_qap_r6();
        assert_eq!(test_circuit.w_mid()[0].evaluate(r6), FE::new(1));
    }

    #[test]
    fn y_mid_test_circuit_on_r5_is_1() {
        let test_circuit = new_test_qap();
        let r5 = Qap::test_qap_r5();
        assert_eq!(test_circuit.y_mid()[0].evaluate(r5), FE::new(1));
    }

    #[test]
    fn y_mid_test_circuit_on_r6_is_0() {
        let test_circuit = new_test_qap();
        let r6 = Qap::test_qap_r6();
        assert_eq!(test_circuit.y_mid()[0].evaluate(r6), FE::new(0));
    }

    #[test]
    fn v_input_test_circuit_has_length_4() {
        let test_circuit = new_test_qap();
        assert_eq!(test_circuit.v_input().len(), 4);
    }

    #[test]
    fn w_input_test_circuit_has_length_4() {
        let test_circuit = new_test_qap();
        assert_eq!(test_circuit.w_input().len(), 4);
    }
    #[test]
    fn y_input_test_circuit_has_length_4() {
        let test_circuit = new_test_qap();
        assert_eq!(test_circuit.y_input().len(), 4);
    }

    #[test]
    fn v_output_test_circuit_has_length_1() {
        let test_circuit = new_test_qap();
        assert_eq!(test_circuit.v_output().len(), 1);
    }

    #[test]
    fn w_output_test_circuit_has_length_1() {
        let test_circuit = new_test_qap();
        assert_eq!(test_circuit.w_output().len(), 1);
    }

    #[test]
    fn y_output_test_circuit_has_length_1() {
        let test_circuit = new_test_qap();
        assert_eq!(test_circuit.y_output().len(), 1);
    }

    #[test]
    /// This test runs multiple cases calculated in paper
    /// t polynomial is tested implicitly by calculating h = p / t
    fn test_polynomial_h_cases() {
        let test_circuit = new_test_qap();

        let inputs = [FE::new(1), FE::new(2), FE::new(3), FE::new(4)];

        let (c5, c6) = Qap::test_qap_solver(inputs);

        let mut c_vector = inputs.to_vec();
        c_vector.append(&mut vec![c5, c6]);

        assert_eq!(
            test_circuit.h_polynomial(&c_vector),
            Polynomial::new_monomial(FE::new(0), 0)
        );

        let inputs = [FE::new(2), FE::new(2), FE::new(2), FE::new(2)];

        let (c5, c6) = Qap::test_qap_solver(inputs);

        let mut c_vector = inputs.to_vec();
        c_vector.append(&mut vec![c5, c6]);

        assert_eq!(
            test_circuit.h_polynomial(&c_vector),
            Polynomial::new_monomial(FE::new(4), 0)
        );

        let inputs = [FE::new(3), FE::new(3), FE::new(3), FE::new(3)];

        let (c5, c6) = Qap::test_qap_solver(inputs);

        let mut c_vector = inputs.to_vec();
        c_vector.append(&mut vec![c5, c6]);

        assert_eq!(
            test_circuit.h_polynomial(&c_vector),
            Polynomial::new_monomial(FE::new(18), 0)
        );

        let inputs = [FE::new(4), FE::new(3), FE::new(2), FE::new(1)];

        let (c5, c6) = Qap::test_qap_solver(inputs);

        let mut c_vector = inputs.to_vec();
        c_vector.append(&mut vec![c5, c6]);

        assert_eq!(
            test_circuit.h_polynomial(&c_vector),
            Polynomial::new_monomial(FE::new(5), 0)
        );
    }

    #[test]
    fn test_circuit_solver_on_2_2_2_2_outputs_4_and_16() {
        let inputs = [FE::new(2), FE::new(2), FE::new(2), FE::new(2)];

        let (c5, c6) = Qap::test_qap_solver(inputs);
        assert_eq!(c5, FE::new(4));
        assert_eq!(c6, FE::new(16));
    }

    #[test]
    fn test_circuit_solver_on_1_2_3_4_outputs_12_and_36() {
        let inputs = [FE::new(1), FE::new(2), FE::new(3), FE::new(4)];

        let (c5, c6) = Qap::test_qap_solver(inputs);
        assert_eq!(c5, FE::new(12));
        assert_eq!(c6, FE::new(36));
    }
}
