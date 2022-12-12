use std::ops;
use rand::prelude::*;

const ORDER: u128 = 18446744073709551359;

#[derive(Debug, PartialEq, Eq)]
pub enum FieldElementError {
    OutOfRangeValue,
    DivisionByZero,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct FieldElement {
    value: u128,
}

impl FieldElement {
    pub fn new(value: u128) -> Result<Self, FieldElementError> {
        if value < ORDER {
            Ok(Self { value })
        } else {
            Err(FieldElementError::OutOfRangeValue)
        }
    }

    pub fn random() -> Self {
        let mut rng = rand::thread_rng();
        let value: u128 = rng.gen();
        // TODO: Some values have higher probability than others.
        FieldElement { value: value % ORDER }
    }

    pub fn pow(self, mut exponent: u128) -> Self {
        let mut result = Self::new(1).unwrap();
        let mut base = self;

        while exponent > 0 {
            // exponent % 2 == 1
            if exponent & 1 == 1 {
                result = result * base;
            }
            // exponent = exponent / 2
            exponent >>= 1;
            base = base * base;
        }
        result
    }

    pub fn inv(self) -> Result<Self, FieldElementError> {
        if self.value != 0 {
            Ok(self.pow(ORDER - 2))
        } else {
            Err(FieldElementError::DivisionByZero)
        }
    }
}

impl ops::Add<FieldElement> for FieldElement {
    type Output = FieldElement;

    fn add(self, a_field_element: FieldElement) -> FieldElement {
        FieldElement::new((self.value + a_field_element.value) % ORDER).unwrap()
    }
}

impl ops::Neg for FieldElement {
    type Output = FieldElement;

    fn neg(self) -> FieldElement {
        FieldElement::new(ORDER - self.value).unwrap()
    }
}

impl ops::Sub<FieldElement> for FieldElement {
    type Output = FieldElement;

    fn sub(self, substrahend: FieldElement) -> FieldElement {
        let neg_substrahend = -substrahend;
        self + neg_substrahend
    }
}

impl ops::Mul for FieldElement {
    type Output = FieldElement;

    fn mul(self, a_field_element: Self) -> Self {
        FieldElement::new(self.value * a_field_element.value % ORDER).unwrap()
    }
}

impl ops::Div for FieldElement {
    type Output = FieldElement;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, dividend: Self) -> Self {
        self * dividend.inv().unwrap()
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn order_must_small_as_to_not_allow_overflows() {
        // ORDER*ORDER < u128::MAX
        assert!(ORDER <= u64::MAX.into());
    }

    #[test]
    fn two_plus_one_is_three() {
        assert_eq!(
            FieldElement::new(2).unwrap() + FieldElement::new(1).unwrap(),
            FieldElement::new(3).unwrap()
        );
    }

    #[test]
    fn max_order_plus_1_is_0() {
        assert_eq!(
            FieldElement::new(ORDER - 1).unwrap() + FieldElement::new(1).unwrap(),
            FieldElement::new(0).unwrap()
        );
    }

    #[test]
    fn when_comparing_13_and_13_they_are_equal() {
        let a: FieldElement = FieldElement::new(13).unwrap();
        let b: FieldElement = FieldElement::new(13).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn when_comparing_13_and_8_they_are_different() {
        let a: FieldElement = FieldElement::new(13).unwrap();
        let b: FieldElement = FieldElement::new(8).unwrap();
        assert_ne!(a, b);
    }

    #[test]
    fn mul_neutral_element() {
        let a: FieldElement = FieldElement::new(1).unwrap();
        let b: FieldElement = FieldElement::new(2).unwrap();
        assert_eq!(a * b, FieldElement::new(2).unwrap());
    }

    #[test]
    fn mul_2_3_is_6() {
        let a: FieldElement = FieldElement::new(2).unwrap();
        let b: FieldElement = FieldElement::new(3).unwrap();
        assert_eq!(a * b, FieldElement::new(6).unwrap());
    }

    #[test]
    fn mul_order_minus_1() {
        let a: FieldElement = FieldElement::new(ORDER - 1).unwrap();
        let b: FieldElement = FieldElement::new(ORDER - 1).unwrap();
        assert_eq!(a * b, FieldElement::new(1).unwrap());
    }

    #[test]
    fn inv_0_error() {
        let a: FieldElement = FieldElement::new(0).unwrap();
        assert_eq!(a.inv().unwrap_err(), FieldElementError::DivisionByZero);
    }

    #[test]
    fn inv_2() {
        let a: FieldElement = FieldElement::new(2).unwrap();
        assert_eq!(a * a.inv().unwrap(), FieldElement::new(1).unwrap());
    }

    #[test]
    fn pow_2_3() {
        assert_eq!(
            FieldElement::new(2).unwrap().pow(3),
            FieldElement::new(8).unwrap()
        )
    }

    #[test]
    fn pow_p_minus_1() {
        assert_eq!(
            FieldElement::new(2).unwrap().pow(ORDER - 1),
            FieldElement::new(1).unwrap()
        )
    }

    #[test]
    fn div_1() {
        assert_eq!(
            FieldElement::new(2).unwrap() / FieldElement::new(1).unwrap(),
            FieldElement::new(2).unwrap()
        )
    }

    #[test]
    fn div_4_2() {
        assert_eq!(
            FieldElement::new(4).unwrap() / FieldElement::new(2).unwrap(),
            FieldElement::new(2).unwrap()
        )
    }

    #[test]
    fn div_4_3() {
        assert_eq!(
            FieldElement::new(4).unwrap() / FieldElement::new(3).unwrap()
                * FieldElement::new(3).unwrap(),
            FieldElement::new(4).unwrap()
        )
    }

    #[test]
    fn two_plus_its_additive_inv_is_0() {
        let two = FieldElement::new(2).unwrap();

        assert_eq!(two + (-two), FieldElement::new(0).unwrap())
    }

    #[test]
    fn four_minus_three_is_1() {
        let four = FieldElement::new(4).unwrap();
        let three = FieldElement::new(3).unwrap();

        assert_eq!(four - three, FieldElement::new(1).unwrap())
    }

    #[test]
    fn zero_minus_1_is_order_minus_1() {
        let zero = FieldElement::new(0).unwrap();
        let one = FieldElement::new(1).unwrap();

        assert_eq!(zero - one, FieldElement::new(ORDER - 1).unwrap())
    }
}
