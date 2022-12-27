mod circuits;
mod math;
mod pinocchio;

use circuits::qap::Qap;
use math::field_element::FieldElement as FE;
use math::polynomial::Polynomial;
use pinocchio::setup::{setup, ToxicWaste};

fn main() {
    let _evaluation_key = setup(Qap::new_test_circuit(), ToxicWaste::sample());

    let element_a = FE::new(8).unwrap();
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

    Qap::new_test_circuit();
    println!("Test cicuit target poly: {:?}", Qap::new_test_circuit());
}
