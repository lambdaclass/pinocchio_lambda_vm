use ark_relations::{
    lc,
    r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError},
};

use crate::fq5::Fq;

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
