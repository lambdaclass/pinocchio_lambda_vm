mod math;
use math::field_element::FieldElement;

fn main() {
    let element_a = FieldElement::new(32);
    let element_b = FieldElement::new(32);
    println!(
        "{:?} + {:?} = {:?}",
        element_a,
        element_b,
        element_a + element_b
    );
}
