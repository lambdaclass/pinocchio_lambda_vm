use super::field_element::FieldElement as FE;
use std::ops;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Polynomial {
    coefficients: Vec<FE>,
}

#[allow(dead_code)] // TODO: Remove this allow.
impl Polynomial {
    pub fn new(coefficients: Vec<FE>) -> Self {
        Polynomial { coefficients }
    }

    fn pad_zeros_like(&mut self, a_polynomial: &Polynomial) {
        let len_difference = a_polynomial.coefficients.len() - self.coefficients.len();
        if len_difference > 0 {
            self.coefficients.extend(vec![FE::zero(); len_difference]);
        }
    }

    fn unpad_zeros(&mut self) {
        while let Some(coefficient) = self.coefficients.last() {
            if *coefficient == FE::zero() {
                self.coefficients.pop();
            } else {
                break;
            }
        }
    }
}

impl ops::Add<Polynomial> for Polynomial {
    type Output = Polynomial;

    fn add(self, a_polynomial: Polynomial) -> Polynomial {
        let coefficients_a = self.coefficients.iter();
        let coefficients_b = a_polynomial.coefficients.iter();
        let new_coefficients = coefficients_a.zip(coefficients_b).map(|(&x, &y)| x + y);
        Polynomial::new(new_coefficients.collect())
    }
}

impl ops::Neg for Polynomial {
    type Output = Polynomial;

    fn neg(self) -> Polynomial {
        Polynomial::new(self.coefficients.iter().map(|&x| -x).collect())
    }
}

impl ops::Sub<Polynomial> for Polynomial {
    type Output = Polynomial;

    fn sub(self, substrahend: Polynomial) -> Polynomial {
        self + (-substrahend)
    }
}

/*
impl ops::Mul for FE {
    type Output = FE;

    fn mul(self, a_field_element: Self) -> Self {
    }
}

impl ops::Div for FE {
    type Output = FE;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, dividend: Self) -> Self {
    }
}
 */

#[cfg(test)]
mod tests {
    /*
        Some of these tests work when the finite field has order greater than 2.
    */
    use super::super::field_element::ORDER;
    use super::*;

    fn polynomial_a() -> Polynomial {
        Polynomial::new(vec![
            FE::new(1).unwrap(),
            FE::new(2).unwrap(),
            FE::new(3).unwrap(),
        ])
    }

    fn polynomial_minus_a() -> Polynomial {
        Polynomial::new(vec![
            FE::new(ORDER - 1).unwrap(),
            FE::new(ORDER - 2).unwrap(),
            FE::new(ORDER - 3).unwrap(),
        ])
    }

    fn polynomial_b() -> Polynomial {
        Polynomial::new(vec![
            FE::new(3).unwrap(),
            FE::new(4).unwrap(),
            FE::new(5).unwrap(),
        ])
    }

    fn polynomial_a_plus_b() -> Polynomial {
        Polynomial::new(vec![
            FE::new(4).unwrap(),
            FE::new(6).unwrap(),
            FE::new(8).unwrap(),
        ])
    }

    fn polynomial_b_minus_a() -> Polynomial {
        Polynomial::new(vec![
            FE::new(2).unwrap(),
            FE::new(2).unwrap(),
            FE::new(2).unwrap(),
        ])
    }

    #[test]
    fn adding_a_and_b_equals_a_plus_b() {
        assert_eq!(polynomial_a() + polynomial_b(), polynomial_a_plus_b());
    }

    #[test]
    fn adding_a_and_a_plus_b_does_not_equal_b() {
        assert_ne!(polynomial_a() + polynomial_a_plus_b(), polynomial_b());
    }

    #[test]
    fn negating_a_is_equal_to_minus_a() {
        assert_eq!(-polynomial_a(), polynomial_minus_a());
    }

    #[test]
    fn negating_a_is_not_equal_to_a() {
        assert_ne!(-polynomial_a(), polynomial_a());
    }

    #[test]
    fn substracting_b_and_a_equals_b_minus_a() {
        assert_eq!(polynomial_b() - polynomial_a(), polynomial_b_minus_a());
    }

    #[test]
    fn unpad_zeros_substract_zeros() {
        let p1 = polynomial_a();
        let mut p2 = polynomial_a();
        p2.coefficients.extend(vec![FE::zero(); 2]);

        assert_ne!(p2, p1);
        p2.unpad_zeros();
        assert_eq!(p2, p1);
    }

    #[test]
    fn pad_zeros_like_adds_zeros_until_matching_size() {
        let mut p1 = polynomial_a();
        let mut p2 = polynomial_a();
        p2.coefficients.extend(vec![FE::zero(); 2]);

        assert_ne!(p2, p1);
        p1.pad_zeros_like(&p2);
        assert_eq!(p2, p1);
    }
}
