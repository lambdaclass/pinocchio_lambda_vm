/// Integration test with the happy path: Setup -> Proof generation -> Proof verification.
use pinocchio_vm::circuits::test_utils::{new_test_r1cs, test_qap_solver};
use pinocchio_vm::math::elliptic_curve::EllipticCurveElement;
use pinocchio_vm::math::field_element::FieldElement as FE;
use pinocchio_vm::pinocchio::prover;
use pinocchio_vm::pinocchio::setup::{setup, EvaluationKey, ToxicWaste, VerificationKey};
use pinocchio_vm::pinocchio::verifier;

fn test_pinocchio(toxic_waste: ToxicWaste) {
    // Get example circuit.
    let test_qap = new_test_r1cs().into();

    // Construct the evaluation and veryfing key.
    let (evaluation_key, verification_key): (
        EvaluationKey<EllipticCurveElement>,
        VerificationKey<EllipticCurveElement>,
    ) = setup(&test_qap, &toxic_waste);

    // Declare inputs to the circuit. Here we choose the ones from the example
    // of the paper.
    let inputs = [FE::new(1), FE::new(2), FE::new(3), FE::new(4)];

    // Execute the circuit with the above inputs.
    // Get the output and all the intermediate values,
    // that is, the values of the circuit wires.
    let (c_mid, c_output) = test_qap_solver(inputs);

    // Construct the witness. In the context of Pinocchio
    // this means the tuple of all the values corresponding
    // to the inputs, the outputs and the values of the
    // circuit wires. The prover needs all this to construct
    // the proof.
    let mut c_vector = inputs.to_vec();
    c_vector.push(c_mid);
    c_vector.push(c_output);

    // Generate the proof of execution.
    let proof = prover::generate_proof(&evaluation_key, &test_qap, &c_vector);

    // Collect all inputs and ouputs of the circuit.
    let mut c_io_vector = inputs.to_vec();
    c_io_vector.push(c_output);

    // Verify the proof.
    let accepted = verifier::verify(&verification_key, &proof, &c_io_vector);

    // Accept or reject the proof.
    assert!(accepted);
}

#[test]
fn test_pinocchio_random_toxic_wate() {
    let toxic_waste = ToxicWaste::sample();
    test_pinocchio(toxic_waste);
}

#[test]
fn test_pinocchio_2() {
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
    test_pinocchio(toxic_waste);
}
