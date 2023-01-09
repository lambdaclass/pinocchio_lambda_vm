use super::{
    cyclic_group::CyclicGroup, field_element::FieldElement,
    field_extension_element::FieldExtensionElement, polynomial::Polynomial,
};

const ORDER_R: u128 = 17; // Base coefficients for polynomials of the circuit (Also the other of the subgroup of the elliptic curve)
const ORDER_P: u128 = 47; // Base coefficients for coordinates of points in elliptic curve
const EMBEDDING_DEGREE: u32 = 4; // Degree to ensure that torsion group is contained in the elliptic curve over field extensions
const ORDER_FIELD_EXTENSION: u128 = ORDER_P.pow(EMBEDDING_DEGREE);
const TARGET_NORMALIZATION_FACTOR: u128 = (ORDER_P.pow(EMBEDDING_DEGREE) - 1) / ORDER_R;
const ELLIPTIC_CURVE_A: u128 = 21;
const ELLIPTIC_CURVE_B: u128 = 15;
const GENERATOR_AFFINE_X: u128 = 45;
const GENERATOR_AFFINE_Y: u128 = 23;

type FEE = FieldExtensionElement;


// Projective Short Weierstrass form
#[derive(Debug, Clone)]
struct EllipticCurveElement {
    x: FEE,
    y: FEE,
    z: FEE,
}

impl EllipticCurveElement {
    fn new(x: FEE, y: FEE, z: FEE) -> Self {
        assert_eq!(Self::defining_equation(&x, &y, &z), FEE::new_base(0));
        Self { x, y, z }
    }

    /// Curve options
    /// Tiny-jubjub curve taken from moonmath (Example 71, page 74) (a=8, b=8, ORDER=13)
    /// Pairing for beginners (Ex. 2.2.2, page 23) (a=905, b=100, ORDER=1021)
    fn defining_equation(x: &FEE, y: &FEE, z: &FEE) -> FEE {
        // println!("{:?} {:?} {:?}", x, y, z);
        y.pow(2) * z - x.pow(3) - Self::a() * x * z.pow(2) - Self::b() * z.pow(3)
    }

    fn a() -> FEE {
        FEE::new_base(ELLIPTIC_CURVE_A)
    }

    fn b() -> FEE {
        FEE::new_base(ELLIPTIC_CURVE_B)
    }

    /// Taken from moonmath (Algorithm 8, page 107)
    fn _pairing(p: &Self, q: &Self) -> FEE {
        if *p == Self::neutral_element() || *q == Self::neutral_element() {
            //|| p == q {
            FEE::new_base(1) // this is equal to [(-1)^r]^(q^k -1)
        } else {
            let (xp, yp) = (&p.x / &p.z, &p.y / &p.z);
            let (xq, yq) = (&q.x / &q.z, &q.y / &q.z);
            let (mut xt, mut yt) = (xp.clone(), yp.clone());
            let (mut f1, mut f2) = (FEE::new_base(1), FEE::new_base(1));

            let mut r = ORDER_R;
            let mut bs = vec![];
            while r > 0 {
                bs.insert(0, r & 1);
                r >>= 1;
            }
            for b in bs[1..].iter() {
                let m = (FEE::new_base(3) * xt.pow(2) + Self::a()) / (FEE::new_base(2) * &yt);

                f1 = &f1.pow(2) * (&yq - &yt - &m * (&xq - &xt));
                f2 = &f2.pow(2) * (&xq + FEE::new_base(2) * &xt - m.pow(2));

                let x2t = m.pow(2) - FEE::new_base(2) * &xt;
                let y2t = -&yt - &m * (&x2t - &xt);
                xt = x2t;
                yt = y2t;

                if *b == 1 {
                    let m = (&yt - &yp) / (&xt - &xp);
                    f1 = f1 * (&yq - &yt - &m * (&xq - &xt));
                    f2 = f2 * (&xq + &xp + &xt - m.pow(2));
                    let xtp = m.pow(2) - &xt - &xp;
                    let ytp = -(&yt) - &m * (&xtp - &xt);
                    xt = xtp;
                    yt = ytp;
                }
            }
            f1 = f1 * (xq - xt);
            (f1 / f2).pow(TARGET_NORMALIZATION_FACTOR)
        }
    }
}

impl PartialEq for EllipticCurveElement {
    fn eq(&self, other: &Self) -> bool {
        // Projective equality relation: first point has to be a multiple of the other
        (&self.x * &other.z == &self.z * &other.x) && (&self.x * &other.y == &self.y * &other.x)
    }
}
impl Eq for EllipticCurveElement {}

impl CyclicGroup for EllipticCurveElement {
    fn generator() -> Self {
        Self::new(
            FEE::new_base(GENERATOR_AFFINE_X),
            FEE::new_base(GENERATOR_AFFINE_Y),
            FEE::new_base(1),
        )
    }

    fn neutral_element() -> Self {
        Self::new(FEE::new_base(0), FEE::new_base(1), FEE::new_base(0))
    }

    fn operate_with_self(&self, times: u128) -> Self {
        let mut times = times;
        let mut result = Self::neutral_element();
        let mut base = self.clone();

        while times > 0 {
            // times % 2 == 1
            if times & 1 == 1 {
                result = result.operate_with(&base);
            }
            // times = times / 2
            times >>= 1;
            base = base.operate_with(&base);
        }
        result
    }

    /// Taken from moonmath (Algorithm 7, page 89)
    fn operate_with(&self, other: &Self) -> Self {
        if *other == Self::neutral_element() {
            self.clone()
        } else if *self == Self::neutral_element() {
            other.clone()
        } else {
            let u1 = &other.y * &self.z;
            let u2 = &self.y * &other.z;
            let v1 = &other.x * &self.z;
            let v2 = &self.x * &other.z;
            if v1 == v2 {
                if u1 != u2 || self.y == FEE::new_base(0) {
                    Self::neutral_element()
                } else {
                    let w = Self::a() * self.z.pow(2) + FEE::new_base(3) * self.x.pow(2);
                    let s = &self.y * &self.z;
                    let b = &self.x * &self.y * &s;
                    let h = w.pow(2) - FEE::new_base(8) * &b;
                    let xp = FEE::new_base(2) * &h * &s;
                    let yp = w * (FEE::new_base(4) * &b - &h)
                        - FEE::new_base(8) * self.y.pow(2) * s.pow(2);
                    let zp = FEE::new_base(8) * s.pow(3);
                    Self::new(xp, yp, zp)
                }
            } else {
                let u = u1 - &u2;
                let v = v1 - &v2;
                let w = &self.z * &other.z;
                let a = u.pow(2) * &w - v.pow(3) - FEE::new_base(2) * v.pow(2) * &v2;
                let xp = &v * &a;
                let yp = u * (v.pow(2) * v2 - a) - v.pow(3) * u2;
                let zp = v.pow(3) * w;
                Self::new(xp, yp, zp)
            }
        }
    }

    fn pairing(&self, _other: &Self) -> Self {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn create_valid_point_works() {
    //     let point = EllipticCurveElement::new(FEE::new_base(13), FEE::new_base(15), FEE::new_base(1));
    //     assert_eq!(point.x, FEE::new_base(13));
    //     assert_eq!(point.y, FEE::new_base(15));
    //     assert_eq!(point.z, FEE::new_base(1));
    //
    //     let point = EllipticCurveElement::new(
    //         FEE::new(Polynomial::new_monomial(FE::new(7), 2)),
    //         FEE::new(Polynomial::new_monomial(FE::new(16), 3)),
    //         FEE::new_base(1)
    //     );
    //     assert_eq!(point.x, FEE::new(Polynomial::new_monomial(FE::new(7), 2)));
    //     assert_eq!(point.y, FEE::new(Polynomial::new_monomial(FE::new(16), 3)));
    //     assert_eq!(point.z, FEE::new_base(1));
    //// }
    //
    //#[test]
    //fn equality_holds_only_for_points_that_are_multiple_of_each_other() {
    //    assert_eq!(
    //        EllipticCurveElement::new(FEE::new_base(26), FEE::new_base(9), FEE::new_base(1)),
    //        EllipticCurveElement::new(FEE::new_base(26), FEE::new_base(9), FEE::new_base(1)),
    //    );
    //
    //    assert_eq!(
    //        EllipticCurveElement::new(FEE::new_base(26), FEE::new_base(9), FEE::new_base(1)),
    //        EllipticCurveElement::new(FEE::new_base(26 * 2), FEE::new_base(9 * 2), FEE::new_base(2)),
    //    );
    //
    //    assert_eq!(
    //        EllipticCurveElement::new(FEE::new_base(0), FEE::new_base(1), FEE::new_base(0)),
    //        EllipticCurveElement::new(FEE::new_base(0), -FEE::new_base(1), FEE::new_base(0)),
    //    );
    //
    //    assert_eq!(
    //        EllipticCurveElement::new(
    //            FEE::new(Polynomial::new_monomial(FE::new(10), 2)),
    //            FEE::new(Polynomial::new_monomial(FE::new(28), 3)),
    //            FEE::new_base(1)),
    //        EllipticCurveElement::new(
    //            FEE::new(Polynomial::new_monomial(FE::new(10), 2)),
    //            FEE::new(Polynomial::new_monomial(FE::new(28), 3)),
    //            FEE::new_base(1)
    //        )
    //    );
    //
    //    assert_ne!(
    //        EllipticCurveElement::new(
    //            FEE::new(Polynomial::new_monomial(FE::new(7), 2)),
    //            FEE::new(Polynomial::new_monomial(FE::new(16), 3)),
    //            FEE::new_base(1)),
    //        EllipticCurveElement::new(
    //            FEE::new(Polynomial::new_monomial(FE::new(10), 2)),
    //            FEE::new(Polynomial::new_monomial(FE::new(28), 3)),
    //            FEE::new_base(1)
    //        )
    //    );
    //
    //    assert_ne!(
    //        EllipticCurveElement::new(FEE::new_base(26), FEE::new_base(9), FEE::new_base(1)),
    //        EllipticCurveElement::new(FEE::new_base(0), FEE::new_base(1), FEE::new_base(0)),
    //    );
    //}
    //
    #[test]
    fn operate_with_self_works() {
        let mut point_1 = EllipticCurveElement::generator();
        point_1 = point_1.operate_with_self(ORDER_R);
        assert_eq!(point_1, EllipticCurveElement::neutral_element());
    }

    //#[test]
    //fn operate_with_self_works_when_using_field_extensions() {
    //    let mut point_1 = EllipticCurveElement::new(
    //        FEE::new(Polynomial::new_monomial(FE::new(17), 2)),
    //        FEE::new(Polynomial::new_monomial(FE::new(15), 3)),
    //        FEE::new_base(1)
    //    );
    //    point_1 = point_1.operate_with_self(ORDER_R);
    //    assert_eq!(point_1, EllipticCurveElement::neutral_element());
    //}
    //
    //#[test]
    //fn test_pairing_between_two_elliptic_curve_elements_over_field_extensions() {
    //    let mut pa = EllipticCurveElement::new(
    //        FEE::new_base(45),
    //        FEE::new_base(23),
    //        FEE::new_base(1)
    //    );
    //    let mut pb = EllipticCurveElement::new(
    //        FEE::new(Polynomial::new(vec![FE::new(29), FE::new(0), FE::new(31)])),
    //        FEE::new(Polynomial::new(vec![FE::new(0), FE::new(11), FE::new(0), FE::new(35)])),
    //        FEE::new_base(1)
    //    );
    //    let mut result = FEE::new(Polynomial::new(vec![
    //        FE::new(39),
    //        FE::new(45),
    //        FE::new(43),
    //        FE::new(33)
    //    ]));
    //
    //    assert_eq!(EllipticCurveElement::_pairing(&pa, &pb), result);
    //}

    #[test]
    #[should_panic]
    fn create_invalid_points_panicks() {
        EllipticCurveElement::new(FEE::new_base(0), FEE::new_base(1), FEE::new_base(1));
    }

    #[test]
    #[should_panic]
    fn create_invalid_points_with_field_extensions_panicks() {
        EllipticCurveElement::new(
            FEE::new(Polynomial::new_monomial(FE::new(1), 2)),
            FEE::new(Polynomial::new_monomial(FE::new(1), 2)),
            FEE::new_base(1),
        );
    }

    #[test]
    fn doubling_a_point_works() {
        println!(
            "{:?}",
            EllipticCurveElement::new(FEE::new_base(45), FEE::new_base(23), FEE::new_base(1))
                .operate_with_self(4)
        )
    }

    #[test]
    fn evaluate_function_on_divisor() {
        let q = EllipticCurveElement::new(
            FEE::new(Polynomial::new(vec![FE::new(29), FE::new(0), FE::new(31)])),
            FEE::new(Polynomial::new(vec![
                FE::new(0),
                FE::new(11),
                FE::new(0),
                FE::new(35),
            ])),
            FEE::new_base(1),
        );
        let double_q = q.operate_with_self(2);
        let double_q = EllipticCurveElement::new(
            double_q.x / &double_q.z,
            double_q.y / &double_q.z,
            FEE::new_base(1),
        );

        let l_div_v_evaluate_double_q =
            (&double_q.y + FEE::new_base(33) * &double_q.x + FEE::new_base(43))
                / &(&double_q.x + FEE::new_base(35));
        let l_div_v_evaluate_q =
            (&q.y + FEE::new_base(33) * &q.x + FEE::new_base(43)) / &(&q.x + FEE::new_base(35));
        let result = &l_div_v_evaluate_double_q / &l_div_v_evaluate_q;

        println!("\n\nq: {:?}", &q);
        println!("2q: {:?}", &double_q);
        println!("l_div_v 2q: {:?}\n\n", &l_div_v_evaluate_double_q);
        println!("l_div_v q: {:?}\n\n", &l_div_v_evaluate_q);
        println!("{:?}\n\n", result);
    }
}
