use super::{
    cyclic_group::CyclicGroup, field_element::FieldElement,
    field_extension_element::FieldExtensionElement, polynomial::Polynomial,
};
use std::ops;


// Curve ("Pairing for beginners", page 57)
const ORDER_R: u128 = 5;
const ORDER_P: u128 = 59;
const EMBEDDING_DEGREE: u32 = 2; // Degree to ensure that torsion group is contained in the elliptic curve over field extensions
const ORDER_FIELD_EXTENSION: u128 = ORDER_P.pow(EMBEDDING_DEGREE);
const TARGET_NORMALIZATION_POWER: u128 = (ORDER_P.pow(EMBEDDING_DEGREE) - 1) / ORDER_R;
const ELLIPTIC_CURVE_A: u128 = 1;
const ELLIPTIC_CURVE_B: u128 = 0;
const GENERATOR_AFFINE_X: u128 = 35;
const GENERATOR_AFFINE_Y: u128 = 31;


type FE = FieldElement<ORDER_P>;
type FEE = FieldExtensionElement;

// Projective Short Weierstrass form
#[derive(Debug, Clone)]
pub struct EllipticCurveElement {
    x: FEE,
    y: FEE,
    z: FEE,
}

impl EllipticCurveElement {
    fn new(x: FEE, y: FEE, z: FEE) -> Self {
        assert_eq!(Self::defining_equation(&x, &y, &z), FEE::new_base(0));
        Self { x, y, z }
    }

    fn defining_equation(x: &FEE, y: &FEE, z: &FEE) -> FEE {
        y.pow(2) * z
            - x.pow(3)
            - FEE::new_base(ELLIPTIC_CURVE_A) * x * z.pow(2)
            - FEE::new_base(ELLIPTIC_CURVE_B) * z.pow(3)
    }

    fn affine(&self) -> Self {
        Self::new(&self.x / &self.z, &self.y / &self.z, FEE::new_base(1))
    }

    fn line(&self, r: &Self, q: &Self) -> FEE {
        assert_ne!(*q, Self::neutral_element());
        if *self == Self::neutral_element() || *r == Self::neutral_element() {
            if self == r {
                return FEE::new_base(1);
            }
            if *self == Self::neutral_element() {
                return &q.x - &r.x;
            } else {
                return &q.x - &self.x;
            }
        } else if self != r {
            if &self.x == &r.x {
                return &q.x - &self.x;
            } else {
                let l = (&r.y - &self.y) / (&r.x - &self.x);
                return &q.y - &self.y - l * (&q.x - &self.x);
            }
        } else {
            let numerator = FEE::new_base(3) * &self.x.pow(2) + FEE::new_base(ELLIPTIC_CURVE_A);
            let denominator = FEE::new_base(2) * &self.y;
            if denominator == FEE::new_base(0) {
                return &q.x - &self.x;
            } else {
                let l = numerator / denominator;
                return &q.y - &self.y - l * (&q.x - &self.x);
            }
        }
    }

    /// Taken from "Pairings for beginners" (Algorithm 5.1, page 79)
    fn miller(p: &Self, q: &Self) -> FEE {
        if *p == Self::neutral_element() || *q == Self::neutral_element() {
            //|| p == q {
            FEE::new_base(1) // this is equal to [(-1)^r]^(q^k -1)
        } else {
            let mut order_r = ORDER_R;
            let mut bs = vec![];
            while order_r > 0 {
                bs.insert(0, order_r & 1);
                order_r >>= 1;
            }

            let mut f = FEE::new_base(1);
            let mut r = p.clone();

            for b in bs[1..].iter() {
                let s = r.operate_with(&r).affine();
                f = f.pow(2) * (r.line(&r, q) / s.line(&-(&s), q));
                r = s;

                if *b == 1 {
                    let mut s = r.operate_with(&p);
                    if s != Self::neutral_element() {
                        s = s.affine();
                    }
                    f = f * (r.line(&p, q) / s.line(&-(&s), q));
                    r = s;
                }
            }
            f
        }
    }

    fn weil_pairing(p: &Self, q: &Self) -> FEE {
        let numerator = Self::miller(p, q);
        let denominator = Self::miller(q, p);
        let result = (numerator / denominator);
        -result
    }

    fn tate_pairing(p: &Self, q: &Self) -> FEE {
        Self::miller(p, q).pow(TARGET_NORMALIZATION_POWER)
    }
    
    fn distorsion_map(p: &Self) -> Self {
        //assert_eq!(Self::defining_equation(&p.x, &p.y, &p.z), FEE::new_base(0));
        //println!("x: {:?} y: {:?} z: {:?}", &p.x, &p.y, &p.z);
        if p == &Self::neutral_element() {
            Self::neutral_element()
        } else {
            let t = FEE::new(Polynomial::new_monomial(FE::new(1), 1));
            let p_affine = p.affine();
            //println!("x: {:?} y: {:?} z: {:?}", -&p.x, &p.y * &t, 1);
            Self::new(-&p_affine.x, &p_affine.y * t, FEE::new_base(1))
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

impl ops::Neg for &EllipticCurveElement {
    type Output = EllipticCurveElement;

    fn neg(self) -> Self::Output {
        Self::Output::new(self.x.clone(), -self.y.clone(), self.z.clone())
    }
}

impl ops::Neg for EllipticCurveElement {
    type Output = EllipticCurveElement;

    fn neg(self) -> Self::Output {
        -&self
    }
}

impl CyclicGroup for EllipticCurveElement {
    type PairingOutput = FEE;

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
                    let w = FEE::new_base(ELLIPTIC_CURVE_A) * self.z.pow(2)
                        + FEE::new_base(3) * self.x.pow(2);
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

    fn pairing(&self, other: &Self) -> Self::PairingOutput {
        let mut p = self.clone();
        let mut q = other.clone();

        if p != Self::neutral_element() {
            p = p.affine();
        }

        if q != Self::neutral_element() {
            q = q.affine();
        }

        Self::tate_pairing(&p, &Self::distorsion_map(&q))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Tiny Jub Jub
    #[test]
    fn create_valid_point_works() {
        let point = EllipticCurveElement::new(FEE::new_base(35), FEE::new_base(31), FEE::new_base(1));
        assert_eq!(point.x, FEE::new_base(35));
        assert_eq!(point.y, FEE::new_base(31));
        assert_eq!(point.z, FEE::new_base(1));
    }

    #[test]
    #[should_panic]
    fn create_invalid_points_panicks() {
        EllipticCurveElement::new(FEE::new_base(0), FEE::new_base(1), FEE::new_base(1));
    }

    #[test]
    fn operate_with_self_works() {
        let mut point_1 = EllipticCurveElement::generator();
        point_1 = point_1.operate_with_self(ORDER_R);
        assert_eq!(point_1, EllipticCurveElement::neutral_element());
    }

    #[test]
    fn doubling_a_point_works() {
        let point = EllipticCurveElement::new(FEE::new_base(35), FEE::new_base(31), FEE::new_base(1));
        let expected_result =
            EllipticCurveElement::new(FEE::new_base(25), FEE::new_base(29), FEE::new_base(1));
        assert_eq!(point.operate_with_self(2).affine(), expected_result);
    }

    #[test]
    fn test_weil_pairing() {
        let mut pa =
            EllipticCurveElement::new(FEE::new_base(35), FEE::new_base(31), FEE::new_base(1));
        let mut pb = EllipticCurveElement::new(
            FEE::new(Polynomial::new(vec![FE::new(24)])),
            FEE::new(Polynomial::new(vec![
                FE::new(0),
                FE::new(31)
            ])),
            FEE::new_base(1),
        );
        let expected_result = FEE::new(Polynomial::new(vec![
            FE::new(46),
            FE::new(3),
        ]));

        let result_weil = EllipticCurveElement::weil_pairing(&pa, &pb);
        assert_eq!(result_weil, expected_result);
    }

    #[test]
    fn test_tate_pairing() {
        let mut pa =
        EllipticCurveElement::new(FEE::new_base(35), FEE::new_base(31), FEE::new_base(1));
        let mut pb = EllipticCurveElement::new(
            FEE::new(Polynomial::new(vec![FE::new(24)])),
            FEE::new(Polynomial::new(vec![
                FE::new(0),
                FE::new(31)
            ])),
            FEE::new_base(1),
        );
        let expected_result = FEE::new(Polynomial::new(vec![
            FE::new(42),
            FE::new(19),
        ]));

        let result_weil = EllipticCurveElement::tate_pairing(&pa, &pb);
        assert_eq!(result_weil, expected_result);
    }
}
