use std::ops;

const ORDER: u64 = 17;

#[derive(Debug)]
pub struct FieldElement {
    value: u64
}


impl FieldElement {
    pub fn new(value: u64) -> Self {
        return Self { value: value % ORDER }
    }
}

impl ops::Add<FieldElement> for FieldElement {
    type Output = FieldElement;

    fn add(self, a_field_element: FieldElement) -> FieldElement {
        return FieldElement::new(self.value + a_field_element.value);
    }
}

impl ops::Mul for FieldElement {
    type Output = FieldElement;

    fn mul(self, a_field_element: Self) -> Self {
        return FieldElement::new(self.value * a_field_element.value);
    }
}

impl Clone for FieldElement {
    fn clone(&self) -> Self {
        return FieldElement::new(self.value);
    }
}

impl PartialEq for FieldElement {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl Eq for FieldElement {}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn adding_13_and_8_when_order_is_17_outputs_4() {
        let a: FieldElement = FieldElement::new(13);
        let b: FieldElement = FieldElement::new(8);
        assert_eq!(a + b, FieldElement::new(4));
    }
}

fn main() {
    
}
