mod bls6_6_fq;

use std::ops::Deref;

use ark_ff::{PrimeField};
use ark_relations::r1cs::ConstraintSystemRef;
use num_bigint::BigUint;
use pinocchio_vm::circuits::r1cs::R1CS;
use pinocchio_vm::math::field_element::FieldElement;
type FE = FieldElement<5>;

pub fn arkworks_cs_to_pinocchio_r1cs<F: PrimeField>(cs: &ConstraintSystemRef<F>) -> R1CS {
    cs.inline_all_lcs();

    let r1cs_matrices = cs.to_matrices().unwrap();

    let a = arkworks_r1cs_matrix_to_pinocchio_r1cs_matrix(
        &r1cs_matrices.a,
        cs.num_witness_variables() + cs.num_instance_variables() - 1,
    );

    let b = arkworks_r1cs_matrix_to_pinocchio_r1cs_matrix(
        &r1cs_matrices.b,
        cs.num_witness_variables() + cs.num_instance_variables() - 1,
    );

    let c = arkworks_r1cs_matrix_to_pinocchio_r1cs_matrix(
        &r1cs_matrices.c,
        cs.num_witness_variables() + cs.num_instance_variables() - 1,
    );

    /* 
        Notice we can't differentiate outputs and inputs from Arkworks CS, but for the proving system everything it just matters it's an IO, so it's a non problem.
    */
    R1CS::new_with_matrixes(a, b, c, cs.num_instance_variables() - 1, 0)
}

pub fn arkworks_witness_and_io_to_pinocchio_witness<F: PrimeField>(cs: &ConstraintSystemRef<F>) -> Vec<FE>{
    let binding = cs.borrow().unwrap();
    let borrowed_cs_ref = binding.deref();
    let arkworks_witness = &borrowed_cs_ref.witness_assignment;
    let arkworks_io = &borrowed_cs_ref.instance_assignment;

    let mut pinocchio_witness = arkworks_io[1..].to_vec();
    
    pinocchio_witness.append(&mut arkworks_witness.to_vec());

    let printable_arkworks_witness: Vec<String> = 
        arkworks_witness.iter().map(|x| x.to_string()).collect();
    println!("Ark witness: {:?}", arkworks_witness);
    println!("Ark witness: {:?}", printable_arkworks_witness);
    //TO DO: Refactor the code to transform a ark Fq to a pinocchio FE
    pinocchio_witness.iter()
        .map(
            |x| 
            biguint_to_u128(x.into_repr().into()))
        .map(|x| (FE::new(x)))
        .collect()
}

fn arkworks_r1cs_matrix_to_pinocchio_r1cs_matrix<F: PrimeField>(
    m: &Vec<Vec<(F, usize)>>,
    total_variables: usize,
) -> Vec<Vec<FE>> {
    sparse_matrix_to_dense(&arkworks_matrix_fps_to_pinocchio_fes(m), total_variables)
}

fn arkworks_matrix_fps_to_pinocchio_fes<F: PrimeField>(
    m: &Vec<Vec<(F, usize)>>,
) -> Vec<Vec<(FE, usize)>> {
    m.iter()
        .map(|x| {
            x.iter()
                .map(|(x, y)| (biguint_to_u128(x.into_repr().into()), *y))
                .map(|(x, y)| (FE::new(x), y))
                .collect()
        })
        .collect()
}

fn sparse_matrix_to_dense(m: &Vec<Vec<(FE, usize)>>, total_variables: usize) -> Vec<Vec<FE>> {
    m.iter()
        .map(|row| sparse_row_to_dense(row, total_variables))
        .collect()
}

fn sparse_row_to_dense(row: &Vec<(FE, usize)>, total_variables: usize) -> Vec<FE> {
    //The first column of the r1cs is used for constants
    // TO DO: Check constants usage

    let mut dense_row = vec![FE::new(0); total_variables + 1];

    // TO DO: Check how constants are set
    dense_row[0] = FE::new(0);
    for element in row {
        println!("Assigned element: {:?}", element.1);
        dense_row[element.1] = element.0;
    }
    dense_row
}

/// Converts biguint to u128
/// If the biguint is bigger than u128 it takes the first 128 bytes
fn biguint_to_u128(big: BigUint) -> u128 {
    match big.to_u64_digits().len() {
        0 => 0,
        1 => big.to_u64_digits()[0].into(),
        _ => big.to_u64_digits()[0] as u128 & ((big.to_u64_digits()[1] as u128) << 64),
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use bls6_6_fq::{Fq};

    use ark_r1cs_std::{fields::fp::FpVar, prelude::AllocVar};
    use ark_relations::{
        lc,
        r1cs::{ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef, SynthesisError},
    };
    use pinocchio_vm::circuits::{r1cs::Constraint, test_utils};
    pub struct MulCircuit {
        /// Public input
        pub a: Fq,
        /// Private input
        pub b: Fq,
    }

    impl ConstraintSynthesizer<Fq> for MulCircuit {
        fn generate_constraints(self, cs: ConstraintSystemRef<Fq>) -> Result<(), SynthesisError> {

            let a = cs.new_witness_variable(|| Ok(self.a))?;

            let b = cs.new_witness_variable(|| Ok(self.b))?;

            let c = cs.new_witness_variable(|| Ok(self.a * self.b))?;
         
            cs.enforce_constraint(lc!() + a, lc!() + b, lc!() + c)?;  
            println!("{:?}*{:?} = {:?}", self.a, self.b, (self.a * self.b).to_string());

            Ok(())
        }
    }

    pub struct TripleMulCircuit {
        /// Public input
        pub a: Fq,
        /// Private input
        pub b: Fq,
        pub c: Fq,
        pub d: Fq,
    }

    impl ConstraintSynthesizer<Fq> for TripleMulCircuit {
        fn generate_constraints(self, cs: ConstraintSystemRef<Fq>) -> Result<(), SynthesisError> {
            let a = FpVar::new_witness(cs.clone(), || Ok(self.a)).unwrap();
            let b = FpVar::new_witness(cs.clone(), || Ok(self.b)).unwrap();
            let c = FpVar::new_witness(cs.clone(), || Ok(self.c)).unwrap();
            let d = FpVar::new_witness(cs.clone(), || Ok(self.d)).unwrap();

            let e = a * b;
            let f = c * d;

            let _g = e * f;

            Ok(())
        }
    }
    pub struct PubMulCircuit {
        /// Public input
        pub a: Fq,
        /// Private input
        pub b: Fq,
    }

    impl ConstraintSynthesizer<Fq> for PubMulCircuit {
        fn generate_constraints(self, cs: ConstraintSystemRef<Fq>) -> Result<(), SynthesisError> {
            let a = cs.new_input_variable(|| Ok(self.a))?;

            let b = cs.new_input_variable(|| Ok(self.b))?;

            let c = cs.new_witness_variable(|| Ok(self.a * self.b))?;

            cs.enforce_constraint(lc!() + a, lc!() + b, lc!() + c)?;

            Ok(())
        }
    }
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

    #[test]
    fn pinocchio_paper_r1cs_from_arkworks_eq_r1cs_from_pinocchio_vm() {
        let a = Fq::from(4);
        let b = Fq::from(3);
        let c = Fq::from(2);
        let d = Fq::from(1);

        let circuit = PinocchioPaperExampleCircuit { a, b, c, d };

        let cs = ConstraintSystem::new_ref();
        circuit.generate_constraints(cs.clone()).unwrap();

        let is_satisfied = cs.is_satisfied().unwrap();
        if !is_satisfied {
            println!("{:?}", cs.which_is_unsatisfied().unwrap().unwrap());
        }

        let converted_r1cs = arkworks_cs_to_pinocchio_r1cs(&cs);

        assert_eq!(
            converted_r1cs.constraints,
            swap_last_variables_test_circuit(&test_utils::new_test_r1cs()).constraints
        );
        assert_eq!(
            converted_r1cs.number_of_inputs + converted_r1cs.number_of_outputs,
            test_utils::new_test_r1cs().number_of_inputs
                + test_utils::new_test_r1cs().number_of_outputs
        );
    }

    #[test]
    fn pinocchio_paper_ark_witness_equals_solver() {
        // 1 + 2 mod 5= 3
        // 3 * 4 mod 5 = 12 mod 5 = 2
        // 3*2 mod 5 = 1
        let a = Fq::from(1);
        //let b = Fq::from((ark_ff::BigInteger256::new([2, 0, 0, 0]));
        let b = Fq::from(2);
        let c = Fq::from(3);
        let d = Fq::from(4);

        let inputs_as_fe = 
        [FE::new(1),FE::new(2),FE::new(3),FE::new(4)];
        let circuit = PinocchioPaperExampleCircuit { a, b, c, d };

        let cs = ConstraintSystem::new_ref();
        circuit.generate_constraints(cs.clone()).unwrap();

        let is_satisfied = cs.is_satisfied().unwrap();

        if !is_satisfied {
            println!("This should be printed");
        }

        
        let ark_witness =         arkworks_witness_and_io_to_pinocchio_witness(&cs);
        
        println!("Witness: {:?}", ark_witness);
        let (c6,c5) = test_utils::test_qap_solver(inputs_as_fe);

        let mut solver_witness = inputs_as_fe.clone().to_vec();
        solver_witness.push(c5);
        solver_witness.push(c6);

        assert_eq!(ark_witness,solver_witness);
    }


    /// This function changes variable 5 for 6
    /// our current implementation of the paper r1cs and the arkworks
    /// version utilizes a different index, but it is the same r1cs
    fn swap_last_variables_test_circuit(r1cs: &R1CS) -> R1CS {
        let mut updated_constraints: Vec<Constraint> = Vec::new();
        for constraint in &r1cs.constraints {
            let mut updated_constraint = constraint.clone();
            updated_constraint.a.swap(5, 6);
            updated_constraint.b.swap(5, 6);
            updated_constraint.c.swap(5, 6);
            updated_constraints.push(updated_constraint);
        }

        R1CS::new(
            updated_constraints,
            r1cs.number_of_inputs,
            r1cs.number_of_outputs,
        )
        .unwrap()
    }

    // The following tests purposes is to see if the transformation
    // is realized without panics
    #[test]
    fn mul_2_3() {
        let a = Fq::from(2);
        let b = Fq::from(3);

        let circuit = MulCircuit { a, b };

        let cs = ConstraintSystem::new_ref();
        circuit.generate_constraints(cs.clone()).unwrap();

        let is_satisfied = cs.is_satisfied().unwrap();
        if !is_satisfied {
            println!("{:?}", cs.which_is_unsatisfied().unwrap().unwrap());
        }
        arkworks_cs_to_pinocchio_r1cs(&cs);
    }

    #[test]
    fn mul_2_3_pub() {
        let a = Fq::from(2);
        let b = Fq::from(3);

        let circuit = PubMulCircuit { a, b };

        let cs = ConstraintSystem::new_ref();
        circuit.generate_constraints(cs.clone()).unwrap();

        let is_satisfied = cs.is_satisfied().unwrap();
        if !is_satisfied {
            println!("{:?}", cs.which_is_unsatisfied().unwrap().unwrap());
        }

        arkworks_cs_to_pinocchio_r1cs(&cs);
    }

    #[test]
    fn triple_mul_1_2_2_2() {
        let a = Fq::from(1);
        let b = Fq::from(2);
        let c = Fq::from(3);
        let d = Fq::from(4);

        let circuit = TripleMulCircuit { a, b, c, d };

        let cs = ConstraintSystem::new_ref();
        circuit.generate_constraints(cs.clone()).unwrap();

        let is_satisfied = cs.is_satisfied().unwrap();
        if !is_satisfied {
            println!("{:?}", cs.which_is_unsatisfied().unwrap().unwrap());
        }

        arkworks_cs_to_pinocchio_r1cs(&cs);
    }
}
