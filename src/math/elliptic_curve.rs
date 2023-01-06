use super::{cyclic_group::CyclicGroup, field_element::FieldElement};

type FE = FieldElement<1021>;
const ORDER: u128 = 23;
// const EMBEDDING_DEGREE: u128 = 11;
const TARGET_NORMALIZATION_FACTOR: u128 = 54645616122141467893743270250140; // (1021 ** EMBEDDING_DEGREE - 1) / ORDER

// Projective Short Weierstrass form
#[derive(Debug, Copy, Clone)]
struct EllipticCurveElement {
    x: FE,
    y: FE,
    z: FE,
}

impl EllipticCurveElement {
    fn new(x: FE, y: FE, z: FE) -> Self {
        assert_eq!(Self::defining_equation(x, y, z), FE::new(0));
        Self { x, y, z }
    }

    /// Curve options
    /// Tiny-jubjub curve taken from moonmath (Example 71, page 74) (a=8, b=8, ORDER=13)
    /// Pairing for beginners (Ex. 2.2.2, page 23) (a=905, b=100, ORDER=1021)
    fn defining_equation(x: FE, y: FE, z: FE) -> FE {
        // println!("{:?} {:?} {:?}", x, y, z);
        y.pow(2) * z - x.pow(3) - Self::a() * x * z.pow(2) - Self::b() * z.pow(3)
    }

    fn a() -> FE {
        FE::new(905)
    }

    fn b() -> FE {
        FE::new(100)
    }

    /// Taken from moonmath (Algorithm 8, page 107)
    fn _pairing(p: Self, q: Self) -> FE {
        if p == Self::neutral_element() || q == Self::neutral_element() || p == q {
            FE::new(1) // this is equal to [(-1)^r]^(q^k -1)
        } else {
            let (xp, yp) = (p.x / p.z, p.y / p.z);
            let (xq, yq) = (q.x / q.z, q.y / q.z);
            let (mut xt, mut yt) = (xp, yp);
            let (mut f1, mut f2) = (FE::new(1), FE::new(1));

            let mut r = ORDER;
            let mut bs = vec![];
            while r > 0 {
                bs.insert(0, r & 1);
                r >>= 1;
            }
            for b in bs[1..].iter() {
                let m = (FE::new(3) * xt.pow(2) + Self::a()) / (FE::new(2) * yt);
                f1 = f1.pow(2) * (yq - yt - m * (xq - xt));
                f2 = f2.pow(2) * (xq + FE::new(2) * xt - m.pow(2));
                let x2t = m.pow(2) - FE::new(2) * xt;
                let y2t = -yt - m * (x2t - xt);
                xt = x2t;
                yt = y2t;
                if *b == 1 {
                    let m = (yt - yp) / (xt - xp);
                    f1 = f1 * (yq - yt - m * (xq - xt));
                    f2 = f2 * (xq + xp + xt - m.pow(2));
                    println!("{:?}", f2);
                    let xtp = m.pow(2) - xt - xp;
                    let ytp = -yt - m * (xtp - xt);
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
        (self.x * other.z == self.z * other.x) && (self.x * other.y == self.y * other.x)
    }
}
impl Eq for EllipticCurveElement {}

impl CyclicGroup for EllipticCurveElement {
    fn generator() -> Self {
        Self::new(FE::new(1006), FE::new(416), FE::new(1))
    }

    fn neutral_element() -> Self {
        Self::new(FE::new(0), FE::new(1), FE::new(0))
    }

    fn operate_with_self(self, times: u128) -> Self {
        let mut times = times;
        let mut result = Self::neutral_element();
        let mut base = self;

        while times > 0 {
            // times % 2 == 1
            if times & 1 == 1 {
                result = result.operate_with(base);
            }
            // times = times / 2
            times >>= 1;
            base = base.operate_with(base);
        }
        result
    }

    /// Taken from moonmath (Algorithm 7, page 89)
    fn operate_with(self, other: Self) -> Self {
        if other == Self::neutral_element() {
            self
        } else if self == Self::neutral_element() {
            other
        } else {
            let u1 = other.y * self.z;
            let u2 = self.y * other.z;
            let v1 = other.x * self.z;
            let v2 = self.x * other.z;
            if v1 == v2 {
                if u1 != u2 || self.y == FE::new(0) {
                    Self::neutral_element()
                } else {
                    let w = Self::a() * self.z.pow(2) + FE::new(3) * self.x.pow(2);
                    let s = self.y * self.z;
                    let b = self.x * self.y * s;
                    let h = w.pow(2) - FE::new(8) * b;
                    let xp = FE::new(2) * h * s;
                    let yp = w * (FE::new(4) * b - h) - FE::new(8) * self.y.pow(2) * s.pow(2);
                    let zp = FE::new(8) * s.pow(3);
                    Self::new(xp, yp, zp)
                }
            } else {
                let u = u1 - u2;
                let v = v1 - v2;
                let w = self.z * other.z;
                let a = u.pow(2) * w - v.pow(3) - FE::new(2) * v.pow(2) * v2;
                let xp = v * a;
                let yp = u * (v.pow(2) * v2 - a) - v.pow(3) * u2;
                let zp = v.pow(3) * w;
                Self::new(xp, yp, zp)
            }
        }
    }

    fn pairing(self, _other: Self) -> Self {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_valid_point_works() {
        let point = EllipticCurveElement::new(FE::new(1006), FE::new(416), FE::new(1));
        assert_eq!(point.x, FE::new(1006));
        assert_eq!(point.y, FE::new(416));
        assert_eq!(point.z, FE::new(1));
    }

    #[test]
    fn equality_holds_only_for_points_that_are_multiple_of_each_other() {
        assert_eq!(
            EllipticCurveElement::new(FE::new(1006), FE::new(416), FE::new(1)),
            EllipticCurveElement::new(FE::new(1006), FE::new(416), FE::new(1)),
        );

        assert_eq!(
            EllipticCurveElement::new(FE::new(1006), FE::new(416), FE::new(1)),
            EllipticCurveElement::new(FE::new(1006 * 2), FE::new(416 * 2), FE::new(2)),
        );

        assert_eq!(
            EllipticCurveElement::new(FE::new(0), FE::new(1), FE::new(0)),
            EllipticCurveElement::new(FE::new(0), -FE::new(1), FE::new(0)),
        );

        assert_ne!(
            EllipticCurveElement::new(FE::new(1006), FE::new(416), FE::new(1)),
            EllipticCurveElement::new(FE::new(0), FE::new(1), FE::new(0)),
        );
    }

    #[test]
    fn operate_with_self_works() {
        let mut point_1 = EllipticCurveElement::generator();
        point_1 = point_1.operate_with_self(23);
        assert_eq!(point_1, EllipticCurveElement::neutral_element());
    }

    #[test]
    #[should_panic]
    fn create_invalid_points_panicks() {
        EllipticCurveElement::new(FE::new(0), FE::new(1), FE::new(1));
    }
}
