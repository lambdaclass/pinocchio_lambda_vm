use super::field_element::FieldElement as FE;
use std::ops;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Polynomial {
    coefficients: Vec<FE>,
}

impl Polynomial {
    /// Creates a new polynomial with the given coefficients
    pub fn new(coefficients: Vec<FE>) -> Self {
        // Removes uneeded 0 coefficients at the end
        let mut unpadded_coefficients = coefficients
            .into_iter()
            .rev()
            .skip_while(|x| *x == FE::zero())
            .collect::<Vec<FE>>();
        unpadded_coefficients.reverse();
        Polynomial {
            coefficients: unpadded_coefficients,
        }
    }

    pub fn zero() -> Self {
        Self::new(Vec::<FE>::new())
    }

    pub fn monomial(coefficient: FE, degree: usize) -> Self {
        let mut coefficients = vec![FE::zero(); degree];
        coefficients.push(coefficient);
        Self::new(coefficients)
    }

    pub fn interpolate(xs: &[FE], ys: &[FE]) -> Polynomial {
        let mut result = Polynomial::zero();

        for (i, y) in ys.iter().enumerate() {
            let mut y_term = Polynomial::new(vec![*y]);
            for (j, x) in xs.iter().enumerate() {
                if i != j {
                    let denominator = Polynomial::new(vec![FE::one() / (xs[i] - *x)]);
                    let numerator = Polynomial::new(vec![-*x, FE::one()]);
                    y_term = &y_term * &(&numerator * &denominator);
                }
            }
            result = result + y_term;
        }
        result
    }

    pub fn evaluate(&self, x: FE) -> FE {
        self.coefficients
            .iter()
            .enumerate()
            .fold(FE::zero(), |acc, (i, &c)| acc + c * x.pow(i as u128))
    }

    pub fn degree(&self) -> usize {
        if !self.coefficients.is_empty() {
            self.coefficients.len() - 1
        } else {
            0
        }
    }

    pub fn last_coefficient(&self) -> FE {
        if let Some(coefficient) = self.coefficients.last() {
            *coefficient
        } else {
            FE::zero()
        }
    }

    /// Returns two new polynomials with the same amount of coefficients
    /// for temporal use
    fn pad_with_zero_coefficients(pa: &Polynomial, pb: &Polynomial) -> (Polynomial, Polynomial) {
        let mut pa = pa.clone();
        let mut pb = pb.clone();

        if pa.coefficients.len() > pb.coefficients.len() {
            pb.coefficients.resize(pa.coefficients.len(), FE::zero());
        } else {
            pa.coefficients.resize(pb.coefficients.len(), FE::zero());
        }
        (pa, pb)
    }
}

impl ops::Add<Polynomial> for Polynomial {
    type Output = Polynomial;

    fn add(self, a_polynomial: Polynomial) -> Polynomial {
        let (pa, pb) = Self::pad_with_zero_coefficients(&self, &a_polynomial);
        let iter_coeff_pa = pa.coefficients.iter();
        let iter_coeff_pb = pb.coefficients.iter();
        let new_coefficients = iter_coeff_pa.zip(iter_coeff_pb).map(|(&x, &y)| x + y);

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

impl ops::Mul<&Polynomial> for &Polynomial {
    type Output = Polynomial;

    fn mul(self, poly: &Polynomial) -> Polynomial {
        let degree = self.degree() + poly.degree();
        let mut coefficients = vec![FE::zero(); degree + 1];

        for i in 0..(poly.degree() + 1) {
            for j in 0..(self.degree() + 1) {
                coefficients[i + j] += poly.coefficients[i] * self.coefficients[j];
            }
        }

        Polynomial::new(coefficients)
    }
}

impl ops::Div<&Polynomial> for &Polynomial {
    type Output = Polynomial;

    fn div(self, dividend: &Polynomial) -> Polynomial {
        if dividend.degree() > self.degree() {
            Polynomial::zero()
        } else {
            let mut n = self.clone();
            let mut q: Vec<FE> = vec![FE::zero(); n.degree() + 1];
            while n != Polynomial::zero() && n.degree() >= dividend.degree() {
                let new_coefficient = n.last_coefficient() / dividend.last_coefficient();
                q[n.degree() - dividend.degree()] = new_coefficient;
                let d = dividend
                    * &Polynomial::monomial(new_coefficient, n.degree() - dividend.degree());
                n = n - d;
            }
            Polynomial::new(q)
        }
    }
}

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
    fn constructor_removes_zeros_at_the_end_of_polynomial() {
        let p1 = Polynomial::new(vec![FE::new(3).unwrap(), FE::new(4).unwrap(), FE::zero()]);
        assert_eq!(
            p1.coefficients,
            vec![FE::new(3).unwrap(), FE::new(4).unwrap()]
        );
    }

    #[test]
    fn pad_with_zero_coefficients_returns_polynomials_with_zeros_until_matching_size() {
        let p1 = Polynomial::new(vec![FE::new(3).unwrap(), FE::new(4).unwrap()]);
        let p2 = Polynomial::new(vec![FE::new(3).unwrap()]);

        assert_eq!(p2.coefficients, vec![FE::new(3).unwrap()]);
        let (pp1, pp2) = Polynomial::pad_with_zero_coefficients(&p1, &p2);
        assert_eq!(pp1, p1);
        assert_eq!(pp2.coefficients, vec![FE::new(3).unwrap(), FE::zero()]);
    }

    #[test]
    fn multiply_2_by_3_is_6() {
        let p1 = Polynomial::new(vec![FE::new(2).unwrap()]);
        let p2 = Polynomial::new(vec![FE::new(3).unwrap()]);
        assert_eq!(&p1 * &p2, Polynomial::new(vec![FE::new(6).unwrap()]));
    }

    #[test]
    fn multiply_2xx_3x_3_times_x_4() {
        let p1 = Polynomial::new(vec![
            FE::new(3).unwrap(),
            FE::new(3).unwrap(),
            FE::new(2).unwrap(),
        ]);
        let p2 = Polynomial::new(vec![FE::new(4).unwrap(), FE::new(1).unwrap()]);
        assert_eq!(
            &p1 * &p2,
            Polynomial::new(vec![
                FE::new(12).unwrap(),
                FE::new(15).unwrap(),
                FE::new(11).unwrap(),
                FE::new(2).unwrap()
            ])
        );
    }

    #[test]
    fn multiply_x_4_times_2xx_3x_3() {
        let p1 = Polynomial::new(vec![
            FE::new(3).unwrap(),
            FE::new(3).unwrap(),
            FE::new(2).unwrap(),
        ]);
        let p2 = Polynomial::new(vec![FE::new(4).unwrap(), FE::new(1).unwrap()]);
        assert_eq!(
            &p2 * &p1,
            Polynomial::new(vec![
                FE::new(12).unwrap(),
                FE::new(15).unwrap(),
                FE::new(11).unwrap(),
                FE::new(2).unwrap()
            ])
        );
    }

    #[test]
    fn division_works() {
        let p1 = Polynomial::new(vec![FE::new(1).unwrap(), FE::new(3).unwrap()]);
        let p2 = Polynomial::new(vec![FE::new(1).unwrap(), FE::new(3).unwrap()]);
        let p3 = &p1 * &p2;
        assert_eq!(&p3 / &p2, p1);
    }

    #[test]
    fn division_by_zero_degree_polynomial_works() {
        let four = FE::new(4).unwrap();
        let two = FE::new(2).unwrap();
        let p1 = Polynomial::new(vec![four, four]);
        let p2 = Polynomial::new(vec![two]);
        assert_eq!(Polynomial::new(vec![two, two]), &p1 / &p2);
    }

    #[test]
    fn evaluate_constant_polynomial_returns_constant() {
        let three = FE::new(3).unwrap();
        let p = Polynomial::new(vec![three]);
        assert_eq!(p.evaluate(FE::new(10).unwrap()), three);
    }

    #[test]
    fn create_degree_0_monomial() {
        assert_eq!(
            Polynomial::monomial(FE::new(3).unwrap(), 0),
            Polynomial::new(vec![FE::new(3).unwrap()])
        );
    }

    #[test]
    fn evaluate_degree_1_monomial() {
        let two = FE::new(2).unwrap();
        let four = FE::new(4).unwrap();
        let p = Polynomial::monomial(two, 1);
        assert_eq!(p.evaluate(two), four);
    }

    #[test]
    fn evaluate_degree_2_monomyal() {
        let two = FE::new(2).unwrap();
        let eight = FE::new(8).unwrap();
        let p = Polynomial::monomial(two, 2);
        assert_eq!(p.evaluate(two), eight);
    }

    #[test]
    fn evaluate_3_term_polynomial() {
        let p = Polynomial::new(vec![
            FE::new(3).unwrap(),
            -FE::new(2).unwrap(),
            FE::new(4).unwrap(),
        ]);
        assert_eq!(p.evaluate(FE::new(2).unwrap()), FE::new(15).unwrap());
    }

    #[test]
    fn simple_interpolating_polynomial_by_hand_works() {
        let denominator = Polynomial::new(vec![
            FE::one() / (FE::new(2).unwrap() - FE::new(4).unwrap()),
        ]);
        let numerator = Polynomial::new(vec![-FE::new(4).unwrap(), FE::one()]);
        let interpolating = &numerator * &denominator;
        assert_eq!(
            (FE::new(2).unwrap() - FE::new(4).unwrap())
                * (FE::one() / (FE::new(2).unwrap() - FE::new(4).unwrap())),
            FE::one()
        );
        assert_eq!(interpolating.evaluate(FE::new(2).unwrap()), FE::one());
        assert_eq!(interpolating.evaluate(FE::new(4).unwrap()), FE::zero());
    }

    #[test]
    fn interpolate_x_2_y_3() {
        let p = Polynomial::interpolate(&[FE::new(2).unwrap()], &[FE::new(3).unwrap()]);
        assert_eq!(FE::new(3).unwrap(), p.evaluate(FE::new(2).unwrap()));
    }

    #[test]
    fn interpolate_x_0_2_y_3_4() {
        let p = Polynomial::interpolate(
            &[FE::new(0).unwrap(), FE::new(2).unwrap()],
            &[FE::new(3).unwrap(), FE::new(4).unwrap()],
        );
        assert_eq!(FE::new(3).unwrap(), p.evaluate(FE::new(0).unwrap()));
        assert_eq!(FE::new(4).unwrap(), p.evaluate(FE::new(2).unwrap()));
    }

    #[test]
    fn interpolate_x_2_5_7_y_10_19_43() {
        let p = Polynomial::interpolate(
            &[
                FE::new(2).unwrap(),
                FE::new(5).unwrap(),
                FE::new(7).unwrap(),
            ],
            &[
                FE::new(10).unwrap(),
                FE::new(19).unwrap(),
                FE::new(43).unwrap(),
            ],
        );

        assert_eq!(FE::new(10).unwrap(), p.evaluate(FE::new(2).unwrap()));
        assert_eq!(FE::new(19).unwrap(), p.evaluate(FE::new(5).unwrap()));
        assert_eq!(FE::new(43).unwrap(), p.evaluate(FE::new(7).unwrap()));
    }
}
