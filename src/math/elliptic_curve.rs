use super::{cyclic_group::CyclicGroup, field_element::FieldElement as FE};

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
        z * y.pow(2) - x.pow(3) - z.pow(2) * Self::a() * x - z.pow(3) * Self::b()
    }

    fn a() -> FE {
        FE::new(905)
    }

    fn b() -> FE {
        FE::new(100)
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

    fn operate_with_self(self, _times: FE) -> Self {
        todo!()
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
    }

    #[test]
    fn operation_on_elliptic_curve_works() {
        let point_1 = EllipticCurveElement::generator();
        let mut res = EllipticCurveElement::neutral_element();

        for _ in 0..23 {
            res = res.operate_with(point_1);
        }
        assert_eq!(res, EllipticCurveElement::neutral_element());
    }

    #[test]
    #[should_panic]
    fn create_invalid_points_panicks() {
        EllipticCurveElement::new(FE::new(0), FE::new(1), FE::new(1));
    }
}
