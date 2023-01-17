use ark_relations::{
    lc,
    r1cs::{ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef, SynthesisError},
};

use pinocchio_vm::{
    circuits::qap::QuadraticArithmeticProgram as Qap,
    math::elliptic_curve::EllipticCurveElement,
    pinocchio::{
        prover,
        setup::{setup, EvaluationKey, ToxicWaste, VerificationKey},
        verifier,
    },
};

use crate::{
    arkworks_cs_to_pinocchio_r1cs, arkworks_io_and_witness_to_pinocchio_io_and_witness, fq5::Fq,
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

    let r1cs = arkworks_cs_to_pinocchio_r1cs(&cs);
    let (io, witness) = arkworks_io_and_witness_to_pinocchio_io_and_witness(&cs);

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

//TO DO: Move this module to a common library
pub struct PinocchioPaperExampleCircuit {
    /// Public input
    pub a: Fq,
    /// Private input
    pub b: Fq,
    pub c: Fq,
    pub d: Fq,
}

impl ConstraintSynthesizer<Fq> for PinocchioPaperExampleCircuit {
    fn generate_constraints(self, cs: ConstraintSystemRef<Fq>) -> Result<(), SynthesisError> {
        let a = cs.new_input_variable(|| Ok(self.a))?;
        let b = cs.new_input_variable(|| Ok(self.b))?;
        let c = cs.new_input_variable(|| Ok(self.c))?;
        let d = cs.new_input_variable(|| Ok(self.d))?;

        let e = cs.new_witness_variable(|| Ok(self.c * self.d))?;
        cs.enforce_constraint(lc!() + c, lc!() + d, lc!() + e)?;

        let calculated_result = self.c * self.d * (self.a + self.b);

        let result = cs.new_input_variable(|| Ok(calculated_result))?;

        cs.enforce_constraint(lc!() + a + b, lc!() + e, lc!() + result)?;

        Ok(())
    }
}
