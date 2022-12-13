mod math;
use math::field_element::FieldElement as FE;
use math::polynomial::interpolate;

fn main() {
    let element_a = FE::new(32).unwrap();
    let element_b = FE::new(32).unwrap();
    println!(
        "{:?} + {:?} = {:?}",
        element_a,
        element_b,
        element_a + element_b
    );

    let x = interpolate(
        &[
            FE::new(0).unwrap(),
            FE::new(1).unwrap(),
            FE::new(2).unwrap(),
        ],
        &[
            FE::new(3).unwrap(),
            FE::new(4).unwrap(),
            FE::new(42).unwrap(),
        ],
        FE::new(3).unwrap(),
    );
    println!("{:?}", x);
}
