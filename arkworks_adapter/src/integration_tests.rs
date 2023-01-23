use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystem};

use pinocchio_lambda_vm::{
    circuits::qap::QuadraticArithmeticProgram as Qap,
    math::elliptic_curve::EllipticCurveElement,
    pinocchio::{
        prover,
        setup::{setup, EvaluationKey, ToxicWaste, VerificationKey},
        verifier,
    },
};

use crate::{
    fq5::Fq, pinocchio_io_and_witness_from_arkworks_cs, pinocchio_r1cs_from_arkworks_cs,
    test_utils::PinocchioPaperExampleCircuit,
};

#[test]
fn create_proof_from_arkworks_and_verify_it() {
    let toxic_waste = ToxicWaste::sample();

    let a = Fq::from(1);
    let b = Fq::from(2);
    let c = Fq::from(3);
    let d = Fq::from(4);

    let circuit = PinocchioPaperExampleCircuit { a, b, c, d };

    let cs = ConstraintSystem::new_ref();
    circuit.generate_constraints(cs.clone()).unwrap();

    let r1cs = pinocchio_r1cs_from_arkworks_cs(&cs);
    let (io, witness) = pinocchio_io_and_witness_from_arkworks_cs(&cs);

    let qap: Qap = r1cs.into();

    let (ek, vk): (
        EvaluationKey<EllipticCurveElement>,
        VerificationKey<EllipticCurveElement>,
    ) = setup(&qap, &toxic_waste);

    let mut c_vector = io.clone();
    c_vector.extend(witness);

    //Reminder: While we differentiate inptus and outputs in Pinocchio
    //All the IO can be placed in the input part, since there is only
    //IO for the proving system
    let proof = prover::generate_proof(&ek, &qap, &c_vector);

    let accepted = verifier::verify(&vk, &proof, &io);

    assert!(accepted);
}
