use super::field_element::FieldElement;

trait GroupElement {
    fn generator() -> Self;
    fn mul_by_scalar(self, other: FieldElement) -> Self;
    fn pairing(self, other: Self) -> Self;
}

impl GroupElement for FieldElement {
    fn generator() -> FieldElement {
        FieldElement::one()
    }

    fn mul_by_scalar(self, other: FieldElement) -> Self {
        self * other
    }

    fn pairing(self, other: Self) -> Self {
        self * other
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn field_element_as_group_element_generator_returns_one() {
        assert_eq!(FieldElement::generator(), FieldElement::one());
    }

    #[test]
    fn field_element_as_group_element_multiplication_by_scalar_works_as_multiplication_in_finite_fields(
    ) {
        let a = FieldElement::new(3).unwrap();
        let b = FieldElement::new(12).unwrap();
        assert_eq!(a * b, a.mul_by_scalar(b));
    }

    #[test]
    fn field_element_as_group_element_pairing_works_as_multiplication_in_finite_fields() {
        let a = FieldElement::new(3).unwrap();
        let b = FieldElement::new(12).unwrap();
        assert_eq!(a * b, a.pairing(b));
    }
}
