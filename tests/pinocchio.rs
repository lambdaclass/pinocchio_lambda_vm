use pinocchio_vm::circuits::qap::{new_test_qap, Qap};
use pinocchio_vm::math::field_element::FieldElement as FE;
use pinocchio_vm::pinocchio::prover;
use pinocchio_vm::pinocchio::setup::{setup, ToxicWaste};
use pinocchio_vm::pinocchio::verifier;

#[test]
fn test_pinocchio() {
    let test_qap = new_test_qap();
    let toxic_waste = ToxicWaste::sample();
    let (evaluation_key, verifying_key) = setup(&test_qap, &toxic_waste);

    let inputs = [
        FE::new(1).unwrap(),
        FE::new(2).unwrap(),
        FE::new(3).unwrap(),
        FE::new(4).unwrap(),
    ];

    // For the test circuit c_mid and c_output has only one element
    let (c_mid, c_output) = Qap::test_qap_solver(inputs);

    let mut c_vector = inputs.to_vec();
    c_vector.push(c_mid);
    c_vector.push(c_output);

    let proof = prover::generate_proof(&evaluation_key, &test_qap, &c_vector);

    let mut c_io_vector = inputs.to_vec();
    c_io_vector.push(c_output);

    let accepted = verifier::verify(&verifying_key, &proof, &c_io_vector);

    assert!(accepted);
}
