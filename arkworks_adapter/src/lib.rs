mod fq5;
#[cfg(test)]
mod integration_tests;
mod test_utils;

use std::ops::Deref;

use ark_ff::PrimeField;
use ark_relations::r1cs::ConstraintSystemRef;
use num_bigint::BigUint;
use pinocchio_lambda_vm::circuits::r1cs::R1CS;
use pinocchio_lambda_vm::config::ORDER_R;

use pinocchio_lambda_vm::math::field_element::FieldElement;
type FE = FieldElement<ORDER_R>;

/// Generates an `R1CS` compatible with Lambda Pinocchio from an Arkworks `ConstraintSystemRef`
/// It supports any ConstraintSystem which isn't using constraints with explicit
/// constants like z=x*y+3
pub fn pinocchio_r1cs_from_arkworks_cs<F: PrimeField>(cs: &ConstraintSystemRef<F>) -> R1CS {
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
        Notice we can't differentiate outputs and inputs from Arkworks CS, but for the proving system everything that matters is that it's public data (IO),
        or private data (witness/c_mid)
    */

    R1CS::new_with_matrixes(a, b, c, cs.num_instance_variables() - 1, 0).unwrap()
}

/// Generates pinocchio IO and Witness from an Arkworks `ConstraintSystemRef`
pub fn pinocchio_io_and_witness_from_arkworks_cs<F: PrimeField>(
    cs: &ConstraintSystemRef<F>,
) -> (Vec<FE>, Vec<FE>) {
    let binding = cs.borrow().unwrap();
    let borrowed_cs_ref = binding.deref();

    let ark_witness = &borrowed_cs_ref.witness_assignment;
    let ark_io = &borrowed_cs_ref.instance_assignment[1..].to_vec();

    let io: Vec<FE> = ark_io.iter().map(ark_fq_to_pinocchio_fe).collect();

    let witness: Vec<FE> = ark_witness.iter().map(ark_fq_to_pinocchio_fe).collect();

    (io, witness)
}

fn arkworks_r1cs_matrix_to_pinocchio_r1cs_matrix<F: PrimeField>(
    m: &[Vec<(F, usize)>],
    total_variables: usize,
) -> Vec<Vec<FE>> {
    sparse_matrix_to_dense(&arkworks_matrix_fps_to_pinocchio_fes(m), total_variables)
}

fn arkworks_matrix_fps_to_pinocchio_fes<F: PrimeField>(
    m: &[Vec<(F, usize)>],
) -> Vec<Vec<(FE, usize)>> {
    m.iter()
        .map(|x| {
            x.iter()
                .map(|(x, y)| (ark_fq_to_pinocchio_fe(x), *y))
                .collect()
        })
        .collect()
}

fn sparse_matrix_to_dense(m: &[Vec<(FE, usize)>], total_variables: usize) -> Vec<Vec<FE>> {
    m.iter()
        .map(|row| sparse_row_to_dense(row, total_variables))
        .collect()
}

fn sparse_row_to_dense(row: &Vec<(FE, usize)>, total_variables: usize) -> Vec<FE> {
    //The first column of the r1cs is used for constants

    let mut dense_row = vec![FE::new(0); total_variables + 1];

    for element in row {
        dense_row[element.1] = element.0;
    }
    dense_row
}

/// Converts an Arkworks fq to a pinocchio FE
fn ark_fq_to_pinocchio_fe<F: PrimeField>(ark_fq: &F) -> FE {
    // into_repr changes back the FQ from the Montgomery to
    // the underlaying representation
    let ark_fq = ark_fq.into_repr();
    let ark_fq_big_int: BigUint = ark_fq.into();
    let fq_as_128 = biguint_to_u128(ark_fq_big_int);
    FE::new(fq_as_128)
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
    use crate::test_utils::PinocchioPaperExampleCircuit;
    use fq5::Fq;

    use ark_relations::{
        lc,
        r1cs::{
            ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef, SynthesisError, Variable,
        },
    };
    use pinocchio_lambda_vm::circuits::{r1cs::Constraint, test_utils};

    // These are tests circuits
    pub struct FullyPrivateMulCircuit {
        pub a: Fq,
        pub b: Fq,
    }

    impl ConstraintSynthesizer<Fq> for FullyPrivateMulCircuit {
        fn generate_constraints(self, cs: ConstraintSystemRef<Fq>) -> Result<(), SynthesisError> {
            let a = cs.new_witness_variable(|| Ok(self.a))?;

            let b = cs.new_witness_variable(|| Ok(self.b))?;

            let c = cs.new_witness_variable(|| Ok(self.a * self.b))?;

            cs.enforce_constraint(lc!() + a, lc!() + b, lc!() + c)?;

            Ok(())
        }
    }
    pub struct PrivInputPubResultMulCircuit {
        pub a: Fq,
        pub b: Fq,
    }

    impl ConstraintSynthesizer<Fq> for PrivInputPubResultMulCircuit {
        fn generate_constraints(self, cs: ConstraintSystemRef<Fq>) -> Result<(), SynthesisError> {
            let a = cs.new_witness_variable(|| Ok(self.a))?;

            let b = cs.new_witness_variable(|| Ok(self.b))?;

            let c = cs.new_input_variable(|| Ok(self.a * self.b))?;

            cs.enforce_constraint(lc!() + a, lc!() + b, lc!() + c)?;

            Ok(())
        }
    }

    // Tests start here
    pub struct MulPlusThree {
        pub a: Fq,
        pub b: Fq,
    }

    impl ConstraintSynthesizer<Fq> for MulPlusThree {
        fn generate_constraints(self, cs: ConstraintSystemRef<Fq>) -> Result<(), SynthesisError> {
            let a = cs.new_witness_variable(|| Ok(self.a))?;

            let b = cs.new_witness_variable(|| Ok(self.b))?;

            let c = cs.new_input_variable(|| Ok(self.a * self.b + Fq::from(3)))?;

            cs.enforce_constraint(
                lc!() + a,
                lc!() + b,
                lc!() + c - (Fq::from(3), Variable::One),
            )?;
            Ok(())
        }
    }

    #[test]
    fn r1cs_from_arkworks_mul_plus_three_first_constraint_c0_is_minus_three() {
        let a = Fq::from(2);
        let b = Fq::from(3);

        let circuit = MulPlusThree { a, b };

        let cs = ConstraintSystem::new_ref();
        circuit.generate_constraints(cs.clone()).unwrap();

        let is_satisfied = cs.is_satisfied().unwrap();
        if !is_satisfied {
            panic!()
        }

        let r1cs = pinocchio_r1cs_from_arkworks_cs(&cs);
        assert_eq!(r1cs.constraints[0].c[0], -FE::new(3));
    }

    #[test]
    fn io_from_arkworks_mul_plus_three_with_0_3_is_3() {
        let a = Fq::from(0);
        let b = Fq::from(3);

        let circuit = MulPlusThree { a, b };

        let cs = ConstraintSystem::new_ref();
        circuit.generate_constraints(cs.clone()).unwrap();

        let (io, _) = pinocchio_io_and_witness_from_arkworks_cs(&cs);

        assert_eq!(io[0], FE::new(3))
    }

    #[test]
    fn io_from_arkworks_mul_plus_three_with_1_2_is_5() {
        let a = Fq::from(1);
        let b = Fq::from(2);

        let circuit = MulPlusThree { a, b };

        let cs = ConstraintSystem::new_ref();
        circuit.generate_constraints(cs.clone()).unwrap();

        let is_satisfied = cs.is_satisfied().unwrap();
        if !is_satisfied {
            panic!()
        }

        let (io, _) = pinocchio_io_and_witness_from_arkworks_cs(&cs);

        assert_eq!(io[0], FE::new(5))
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
            panic!()
        }

        let converted_r1cs = pinocchio_r1cs_from_arkworks_cs(&cs);

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
    fn pinocchio_paper_ark_witness_and_ark_io_are_correct() {
        let a = Fq::from(1);
        let b = Fq::from(2);
        let c = Fq::from(3);
        let d = Fq::from(4);

        let inputs_as_fe = [FE::new(1), FE::new(2), FE::new(3), FE::new(4)];
        let circuit = PinocchioPaperExampleCircuit { a, b, c, d };

        let cs = ConstraintSystem::new_ref();
        circuit.generate_constraints(cs.clone()).unwrap();

        let is_satisfied = cs.is_satisfied().unwrap();

        if !is_satisfied {
            panic!()
        }

        let (ark_io, ark_witness) = pinocchio_io_and_witness_from_arkworks_cs(&cs);

        //Since the Arcworks version of the circuit has the last variables
        //switched, c6 and c5 needs to be switched
        let (c5, c6) = test_utils::test_qap_solver(inputs_as_fe);

        let solver_witness = vec![c5];
        let mut io = inputs_as_fe.to_vec();
        io.push(c6);

        assert_eq!(ark_witness, solver_witness);
        assert_eq!(ark_io, io);
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

    #[test]
    fn fully_private_mul_circuit_has_1_constraint_in_pinocchio_r1cs() {
        let a = Fq::from(2);
        let b = Fq::from(3);

        let circuit = FullyPrivateMulCircuit { a, b };

        let cs = ConstraintSystem::new_ref();
        circuit.generate_constraints(cs.clone()).unwrap();

        let is_satisfied = cs.is_satisfied().unwrap();
        if !is_satisfied {
            panic!()
        }
        let r1cs = pinocchio_r1cs_from_arkworks_cs(&cs);
        assert_eq!(r1cs.number_of_constraints(), 1)
    }

    #[test]
    fn fully_private_mul_circuit_has_3_elements_in_witness_and_0_io_for_pinocchio() {
        let a = Fq::from(2);
        let b = Fq::from(3);

        let circuit = FullyPrivateMulCircuit { a, b };

        let cs = ConstraintSystem::new_ref();
        circuit.generate_constraints(cs.clone()).unwrap();

        let is_satisfied = cs.is_satisfied().unwrap();
        if !is_satisfied {
            panic!()
        }
        let (io, witness) = pinocchio_io_and_witness_from_arkworks_cs(&cs);
        assert_eq!(io.len(), 0);
        assert_eq!(witness.len(), 3);
    }

    #[test]
    fn priv_input_pub_output_mul_has_1_io_and_2_witness() {
        let a = Fq::from(2);
        let b = Fq::from(3);

        let circuit = PrivInputPubResultMulCircuit { a, b };

        let cs = ConstraintSystem::new_ref();
        circuit.generate_constraints(cs.clone()).unwrap();

        let is_satisfied = cs.is_satisfied().unwrap();
        if !is_satisfied {
            panic!()
        }

        let (io, witness) = pinocchio_io_and_witness_from_arkworks_cs(&cs);

        assert_eq!(io.len(), 1);
        assert_eq!(witness.len(), 2);
    }

    #[test]
    fn priv_input_pub_output_mul_2_3_in_pinocchio_has_io_eq_to_6() {
        let a = Fq::from(2);
        let b = Fq::from(3);

        let circuit = PrivInputPubResultMulCircuit { a, b };

        let cs = ConstraintSystem::new_ref();
        circuit.generate_constraints(cs.clone()).unwrap();

        let is_satisfied = cs.is_satisfied().unwrap();
        if !is_satisfied {
            panic!()
        }

        let (io, _) = pinocchio_io_and_witness_from_arkworks_cs(&cs);

        assert_eq!(io, [FE::new(6)]);
    }
}
