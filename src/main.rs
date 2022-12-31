mod circuits;
mod math;
mod pinocchio;

use circuits::qap::Qap;
use pinocchio::prover;
use pinocchio::setup::{setup, ToxicWaste};

fn main() {
    let test_qap = Qap::new_test_qap();
    let toxic_waste = ToxicWaste::sample();
    let (evaluation_key, verifying_key) = setup(&test_qap, &toxic_waste);
    let proof = prover::generate_proof(&evaluation_key, &test_qap);
    println!("Evaluation key: {:?}", evaluation_key);
    println!("Verifying key: {:?}", verifying_key);
    println!("Proof: {:?}", proof);
}
