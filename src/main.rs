mod circuits;
mod math;
mod pinocchio;

use circuits::qap::Qap;
use math::field_element::FieldElement as FE;
use pinocchio::prover;
use pinocchio::setup::{setup, ToxicWaste};

fn main() {
    let test_qap = Qap::new_test_qap();
    let toxic_waste = ToxicWaste::sample();
    let (evaluation_key, verifying_key) = setup(&test_qap, &toxic_waste);

    let inputs = [
        FE::new(3).unwrap(),
        FE::new(3).unwrap(),
        FE::new(3).unwrap(),
        FE::new(3).unwrap(),
    ];

    // For the test circuit c_mid and c_output has only one element
    let (c_mid, c_output) = Qap::test_qap_solver(inputs);

    let mut c_vector = vec![FE::zero()];
    c_vector.append(&mut inputs.to_vec());
    c_vector.push(c_mid);
    c_vector.push(c_output);

    let proof = prover::generate_proof(&evaluation_key, &test_qap, &c_vector);

    println!("Evaluation key: {:?}", evaluation_key);
    println!("Verifying key: {:?}", verifying_key);
    println!("Proof: {:?}", proof);
}
