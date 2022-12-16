mod math;
use math::field_element::FieldElement as FE;
use math::polynomial::Polynomial;

fn main() {
    let element_a = FE::new(32).unwrap();
    let element_b = FE::new(32).unwrap();
    println!(
        "{:?} + {:?} = {:?}",
        element_a,
        element_b,
        element_a + element_b
    );

    let y = Polynomial::interpolate(
        &[
            FE::new(1).unwrap(),
            FE::new(2).unwrap(),
            FE::new(3).unwrap(),
        ],
        &[
            FE::new(1).unwrap(),
            FE::new(4).unwrap(),
            FE::new(8).unwrap(),
        ],
    );
    println!("poly = {:?}", y);
    println!("poly(3) = {:?}", y.evaluate(FE::new(3).unwrap()));
}
