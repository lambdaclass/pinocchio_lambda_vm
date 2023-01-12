use super::{field_element::FieldElement, polynomial::Polynomial};
use std::ops;

const ORDER_P: u128 = 59; // Base coefficients for coordinates of points in elliptic curve
const EMBEDDING_DEGREE: u32 = 2; // Degree to ensure that torsion group is contained in the elliptic curve over field extensions
const ORDER_FIELD_EXTENSION: u128 = ORDER_P.pow(EMBEDDING_DEGREE);
type FE = FieldElement<ORDER_P>;

#[derive(Debug, Clone)]
pub struct FieldExtensionElement {
    value: Polynomial<ORDER_P>,
}

// Taken from Moonmath BLS6_6 curve (Page 129)
impl FieldExtensionElement {
    pub fn new(value: Polynomial<ORDER_P>) -> Self {
        let (_quotient, remainder) =
            value.long_division_with_remainder(&Self::defining_polynomial());
        Self { value: remainder }
    }

    pub fn new_base(value: u128) -> Self {
        Self::new(Polynomial::new_monomial(FE::new(value), 0))
    }

    pub fn defining_polynomial() -> Polynomial<ORDER_P> {
        let linear_term = Polynomial::new_monomial(FE::new(1), 0);
        let higher_order_term = Polynomial::new_monomial(FE::new(1), 2);
        linear_term + higher_order_term
    }

    pub fn pow(&self, mut exponent: u128) -> Self {
        let mut result = Self::new(Polynomial::new_monomial(FE::new(1), 0));
        let mut base = self.clone();

        while exponent > 0 {
            // exponent % 2 == 1
            if exponent & 1 == 1 {
                result = &result * &base;
            }
            // exponent = exponent / 2
            exponent >>= 1;
            base = &base * &base;
        }
        result
    }

    pub fn add(&self, other: &Self) -> Self {
        self + other
    }

    pub fn sub(&self, other: &Self) -> Self {
        self - other
    }

    pub fn mul(&self, other: &Self) -> Self {
        self * other
    }

    pub fn div(&self, other: &Self) -> Self {
        self / other
    }

    pub fn inv(self) -> Self {
        assert_ne!(self.value, Polynomial::zero());
        self.pow(ORDER_FIELD_EXTENSION - 2)
    }
}

impl PartialEq for FieldExtensionElement {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}
impl Eq for FieldExtensionElement {}

impl ops::Add<&FieldExtensionElement> for &FieldExtensionElement {
    type Output = FieldExtensionElement;

    fn add(self, a_field_element: &FieldExtensionElement) -> Self::Output {
        Self::Output {
            value: &a_field_element.value + &self.value,
        }
    }
}

impl ops::Add<FieldExtensionElement> for FieldExtensionElement {
    type Output = FieldExtensionElement;

    fn add(self, a_field_element: FieldExtensionElement) -> Self::Output {
        &self + &a_field_element
    }
}

impl ops::Add<&FieldExtensionElement> for FieldExtensionElement {
    type Output = FieldExtensionElement;

    fn add(self, a_field_element: &FieldExtensionElement) -> Self::Output {
        &self + a_field_element
    }
}

impl ops::Add<FieldExtensionElement> for &FieldExtensionElement {
    type Output = FieldExtensionElement;

    fn add(self, a_field_element: FieldExtensionElement) -> Self::Output {
        self + &a_field_element
    }
}

impl ops::Neg for &FieldExtensionElement {
    type Output = FieldExtensionElement;

    fn neg(self) -> Self::Output {
        Self::Output {
            value: -self.value.clone(),
        }
    }
}

impl ops::Neg for FieldExtensionElement {
    type Output = FieldExtensionElement;

    fn neg(self) -> Self::Output {
        -&self
    }
}

impl ops::Sub<&FieldExtensionElement> for &FieldExtensionElement {
    type Output = FieldExtensionElement;

    fn sub(self, substrahend: &FieldExtensionElement) -> Self::Output {
        self + &(-substrahend)
    }
}

impl ops::Sub<FieldExtensionElement> for FieldExtensionElement {
    type Output = FieldExtensionElement;

    fn sub(self, substrahend: FieldExtensionElement) -> Self::Output {
        &self - &substrahend
    }
}

impl ops::Sub<&FieldExtensionElement> for FieldExtensionElement {
    type Output = FieldExtensionElement;

    fn sub(self, substrahend: &FieldExtensionElement) -> Self::Output {
        &self - substrahend
    }
}

impl ops::Sub<FieldExtensionElement> for &FieldExtensionElement {
    type Output = FieldExtensionElement;

    fn sub(self, substrahend: FieldExtensionElement) -> Self::Output {
        self - &substrahend
    }
}

impl ops::Mul<&FieldExtensionElement> for &FieldExtensionElement {
    type Output = FieldExtensionElement;

    fn mul(self, a_field_extension_element: &FieldExtensionElement) -> Self::Output {
        let p = self.value.mul_with_ref(&a_field_extension_element.value);
        Self::Output::new(p)
    }
}

impl ops::Mul<FieldExtensionElement> for FieldExtensionElement {
    type Output = FieldExtensionElement;

    fn mul(self, a_field_extension_element: FieldExtensionElement) -> Self::Output {
        &self * &a_field_extension_element
    }
}

impl ops::Mul<&FieldExtensionElement> for FieldExtensionElement {
    type Output = FieldExtensionElement;

    fn mul(self, a_field_extension_element: &FieldExtensionElement) -> Self::Output {
        &self * a_field_extension_element
    }
}

impl ops::Mul<FieldExtensionElement> for &FieldExtensionElement {
    type Output = FieldExtensionElement;

    fn mul(self, a_field_extension_element: FieldExtensionElement) -> Self::Output {
        self * &a_field_extension_element
    }
}

impl ops::Div<&FieldExtensionElement> for &FieldExtensionElement {
    type Output = FieldExtensionElement;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, dividend: &FieldExtensionElement) -> Self::Output {
        self * &dividend.clone().inv()
    }
}

impl ops::Div<FieldExtensionElement> for FieldExtensionElement {
    type Output = FieldExtensionElement;

    fn div(self, dividend: FieldExtensionElement) -> Self::Output {
        &self / &dividend
    }
}

impl ops::Div<FieldExtensionElement> for &FieldExtensionElement {
    type Output = FieldExtensionElement;

    fn div(self, dividend: FieldExtensionElement) -> Self::Output {
        self / &dividend
    }
}

impl ops::Div<&FieldExtensionElement> for FieldExtensionElement {
    type Output = FieldExtensionElement;

    fn div(self, dividend: &FieldExtensionElement) -> Self::Output {
        &self / dividend
    }
}
