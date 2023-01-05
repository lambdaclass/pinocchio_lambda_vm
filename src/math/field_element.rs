use super::cyclic_group::CyclicGroup;
use rand::prelude::*;
use std::ops;

pub const ORDER: u128 = 18446744073709551359;

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
    /// Creates a new field element with `value` modulo order of the field
    pub fn new(value: u128) -> Self {
        Self {
            value: value % ORDER,
        }
    }

    pub fn random() -> Self {
        let value: u128 = rand::thread_rng().gen_range(0..ORDER);
        FieldElement { value }
    }

    pub fn pow(self, mut exponent: u128) -> Self {
        let mut result = Self::new(1);
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
        FieldElement::new(self.value + a_field_element.value)
    }
}

impl ops::AddAssign for FieldElement {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl ops::Neg for FieldElement {
    type Output = FieldElement;

    fn neg(self) -> FieldElement {
        FieldElement::new(ORDER - self.value)
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
        FieldElement::new(self.value * a_field_element.value)
    }
}

impl ops::Div for FieldElement {
    type Output = FieldElement;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, dividend: Self) -> Self {
        self * dividend.inv().unwrap()
    }
}

impl CyclicGroup for FieldElement {
    fn generator() -> FieldElement {
        FieldElement::new(1)
    }

    fn neutral_element() -> FieldElement {
        FieldElement::new(0)
    }

    fn operate_with_self(self, other: FieldElement) -> Self {
        self * other
    }

    fn pairing(self, other: Self) -> Self {
        self * other
    }

    fn operate_with(self, other: Self) -> Self {
        self + other
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
            FieldElement::new(2) + FieldElement::new(1),
            FieldElement::new(3)
        );
    }

    #[test]
    fn max_order_plus_1_is_0() {
        assert_eq!(
            FieldElement::new(ORDER - 1) + FieldElement::new(1),
            FieldElement::new(0)
        );
    }

    #[test]
    fn when_comparing_13_and_13_they_are_equal() {
        let a: FieldElement = FieldElement::new(13);
        let b: FieldElement = FieldElement::new(13);
        assert_eq!(a, b);
    }

    #[test]
    fn when_comparing_13_and_8_they_are_different() {
        let a: FieldElement = FieldElement::new(13);
        let b: FieldElement = FieldElement::new(8);
        assert_ne!(a, b);
    }

    #[test]
    fn mul_neutral_element() {
        let a: FieldElement = FieldElement::new(1);
        let b: FieldElement = FieldElement::new(2);
        assert_eq!(a * b, FieldElement::new(2));
    }

    #[test]
    fn mul_2_3_is_6() {
        let a: FieldElement = FieldElement::new(2);
        let b: FieldElement = FieldElement::new(3);
        assert_eq!(a * b, FieldElement::new(6));
    }

    #[test]
    fn mul_order_minus_1() {
        let a: FieldElement = FieldElement::new(ORDER - 1);
        let b: FieldElement = FieldElement::new(ORDER - 1);
        assert_eq!(a * b, FieldElement::new(1));
    }

    #[test]
    fn inv_0_error() {
        let a: FieldElement = FieldElement::new(0);
        assert_eq!(a.inv().unwrap_err(), FieldElementError::DivisionByZero);
    }

    #[test]
    fn inv_2() {
        let a: FieldElement = FieldElement::new(2);
        assert_eq!(a * a.inv().unwrap(), FieldElement::new(1));
    }

    #[test]
    fn pow_2_3() {
        assert_eq!(FieldElement::new(2).pow(3), FieldElement::new(8))
    }

    #[test]
    fn pow_p_minus_1() {
        assert_eq!(FieldElement::new(2).pow(ORDER - 1), FieldElement::new(1))
    }

    #[test]
    fn div_1() {
        assert_eq!(
            FieldElement::new(2) / FieldElement::new(1),
            FieldElement::new(2)
        )
    }

    #[test]
    fn div_4_2() {
        assert_eq!(
            FieldElement::new(4) / FieldElement::new(2),
            FieldElement::new(2)
        )
    }

    #[test]
    fn div_4_3() {
        assert_eq!(
            FieldElement::new(4) / FieldElement::new(3) * FieldElement::new(3),
            FieldElement::new(4)
        )
    }

    #[test]
    fn two_plus_its_additive_inv_is_0() {
        let two = FieldElement::new(2);

        assert_eq!(two + (-two), FieldElement::new(0))
    }

    #[test]
    fn four_minus_three_is_1() {
        let four = FieldElement::new(4);
        let three = FieldElement::new(3);

        assert_eq!(four - three, FieldElement::new(1))
    }

    #[test]
    fn zero_minus_1_is_order_minus_1() {
        let zero = FieldElement::new(0);
        let one = FieldElement::new(1);

        assert_eq!(zero - one, FieldElement::new(ORDER - 1))
    }

    #[test]
    fn neg_zero_is_zero() {
        let zero = FieldElement::new(0);

        assert_eq!(-zero, zero);
    }

    #[test]
    fn zero_constructor_returns_zero() {
        assert_eq!(FieldElement::new(0), FieldElement::new(0));
    }

    #[test]
    fn field_element_as_group_element_generator_returns_one() {
        assert_eq!(FieldElement::generator(), FieldElement::new(1));
    }

    #[test]
    fn field_element_as_group_element_multiplication_by_scalar_works_as_multiplication_in_finite_fields(
    ) {
        let a = FieldElement::new(3);
        let b = FieldElement::new(12);
        assert_eq!(a * b, a.operate_with_self(b));
    }

    #[test]
    fn field_element_as_group_element_pairing_works_as_multiplication_in_finite_fields() {
        let a = FieldElement::new(3);
        let b = FieldElement::new(12);
        assert_eq!(a * b, a.pairing(b));
    }
}
