use std::ops;

const ORDER: u128 = 1000000;

#[derive(Debug)]
pub enum FieldElementError {
    OutOfRangeValue,
}

#[derive(Debug, Copy, Clone)]
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
}

impl ops::Add<FieldElement> for FieldElement {
    type Output = FieldElement;

    fn add(self, a_field_element: FieldElement) -> FieldElement {
        FieldElement::new(self.value + a_field_element.value).unwrap()
    }
}

impl ops::Mul for FieldElement {
    type Output = FieldElement;

    fn mul(self, a_field_element: Self) -> Self {
        FieldElement::new(self.value * a_field_element.value).unwrap()
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
    fn order_must_small_as_to_not_allow_overflows() {
        // ORDER*ORDER < u128::MAX
        assert!(ORDER <= u64::MAX.into());
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
}
