use pinocchio_vm::circuits::test_utils::{new_test_qap, new_test_r1cs, test_qap_solver};
use pinocchio_vm::math::elliptic_curve::EllipticCurveElement;
use pinocchio_vm::math::field_element::FieldElement as FE;
use pinocchio_vm::pinocchio::prover;
use pinocchio_vm::pinocchio::setup::{setup, EvaluationKey, ToxicWaste, VerifyingKey};
use pinocchio_vm::pinocchio::verifier;

#[test]
fn test_pinocchio() {
    let test_qap = new_test_r1cs().into();
    let toxic_waste = ToxicWaste::sample();
    let (ek, verifying_key): (
        EvaluationKey<EllipticCurveElement>,
        VerifyingKey<EllipticCurveElement>,
    ) = setup(&test_qap, &toxic_waste);

    let inputs = [FE::new(1), FE::new(2), FE::new(3), FE::new(4)];

    // For the test circuit c_mid and c_output has only one element
    let (c_mid, c_output) = test_qap_solver(inputs);

    let mut c_vector = inputs.to_vec();
    c_vector.push(c_mid);
    c_vector.push(c_output);

    let proof = prover::generate_proof(&ek, &test_qap, &c_vector);

    let mut c_io_vector = inputs.to_vec();
    c_io_vector.push(c_output);

    let accepted = verifier::verify(&verifying_key, &proof, &c_io_vector);

    assert!(accepted);
}

#[test]
fn test_pinocchio_2() {
    let test_qap = new_test_qap();
    let toxic_waste = ToxicWaste::new(
        FE::new(2),
        FE::new(3),
        FE::new(3),
        FE::new(4),
        FE::new(2),
        FE::new(3),
        FE::new(2),
        FE::new(3),
    );
    let (ek, vk): (
        EvaluationKey<EllipticCurveElement>,
        VerifyingKey<EllipticCurveElement>,
    ) = setup(&test_qap, &toxic_waste);

    let inputs = [FE::new(1), FE::new(2), FE::new(3), FE::new(4)];

    // For the test circuit c_mid and c_output has only one element
    let (c_mid, c_output) = test_qap_solver(inputs);

    let mut c_vector = inputs.to_vec();
    c_vector.push(c_mid);
    c_vector.push(c_output);

    let proof = prover::generate_proof(&ek, &test_qap, &c_vector);

    let mut c_io_vector = inputs.to_vec();
    c_io_vector.push(c_output);

    let accepted = verifier::verify(&vk, &proof, &c_io_vector);

    assert!(accepted);
}
