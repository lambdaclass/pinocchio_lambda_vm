use super::field_element::FieldElement;

pub trait Group {
    fn generator() -> Self;
    fn mul_by_scalar(self, other: FieldElement) -> Self;
    fn pairing(self, other: Self) -> Self;
}
