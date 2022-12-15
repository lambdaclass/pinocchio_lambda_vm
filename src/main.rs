mod math;
use math::field_element::FieldElement as FE;

fn main() {
    let element_a = FE::new(32).unwrap();
    let element_b = FE::new(32).unwrap();
    let c = FE::zero();
    println!(
        "{:?} + {:?} = {:?}",
        element_a,
        element_b,
        element_a + element_b + c
    );
}
