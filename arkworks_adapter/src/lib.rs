use ark_bn254::{Fq, FqParameters};
use ark_ff::{BigInteger256, Fp256, FpParameters, PrimeField};
use ark_relations::r1cs::ConstraintSystemRef;
use num_bigint::BigUint;
use pinocchio_vm::circuits::r1cs::R1CS;
use pinocchio_vm::math::field_element::FieldElement;

type FE = FieldElement<5>;

pub fn arcworks_cs_to_pinocchio_r1cs(cs: &ConstraintSystemRef<Fp256<FqParameters>>) {
    cs.inline_all_lcs();

    let r1cs_matrices = cs.to_matrices().unwrap();

    println!("CS Wit vars: {:?}", cs.num_witness_variables());
    println!("CS Inst vars: {:?}", cs.num_instance_variables());

    let A = arcworks_r1cs_matrix_to_pinocchio_r1cs_matrix(
        &r1cs_matrices.a,
        cs.num_instance_variables(),
    );

    let B = arcworks_r1cs_matrix_to_pinocchio_r1cs_matrix(
        &r1cs_matrices.b,
        cs.num_instance_variables(),
    );

    let C = arcworks_r1cs_matrix_to_pinocchio_r1cs_matrix(
        &r1cs_matrices.c,
        cs.num_instance_variables(),
    );
}

fn arcworks_r1cs_matrix_to_pinocchio_r1cs_matrix(
    m: &Vec<Vec<(Fp256<FqParameters>, usize)>>,
    num_instance_vars: usize,
) -> Vec<Vec<FE>> {
    sparse_matrix_to_dense(&arcworks_matrix_fps_to_pinocchio_fes(m), num_instance_vars)
}

fn arcworks_matrix_fps_to_pinocchio_fes(
    m: &Vec<Vec<(Fp256<FqParameters>, usize)>>,
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

fn sparse_matrix_to_dense(m: &Vec<Vec<(FE, usize)>>, num_ins_vars: usize) -> Vec<Vec<FE>> {
    m.iter()
        .map(|row| sparse_row_to_dense(row, num_ins_vars))
        .collect()
}

fn sparse_row_to_dense(row: &Vec<(FE, usize)>, num_instance_vars: usize) -> Vec<FE> {
    //The first column of the r1cs is used for constants
    // TO DO: Check constants usage
    let mut dense_row = vec![FE::new(0); num_instance_vars + 1];
    for element in row {
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
    use ark_bn254::Fq;
    use ark_ff::BigInteger256;
    use ark_relations::{
        lc,
        r1cs::{ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef, SynthesisError},
    };

    use super::*;

    pub struct MulCircuit {
        /// Public input
        pub a: Fq,
        /// Private input
        pub b: Fq,
    }

    impl ConstraintSynthesizer<Fq> for MulCircuit {
        fn generate_constraints(self, cs: ConstraintSystemRef<Fq>) -> Result<(), SynthesisError> {
            let a = cs.new_input_variable(|| Ok(self.a))?;

            let b = cs.new_input_variable(|| Ok(self.b))?;

            let c = cs.new_witness_variable(|| Ok(self.b))?;

            //c = a.mul(b);
            cs.enforce_constraint(lc!() + a, lc!() + b, lc!() + c)?;

            Ok(())
        }
    }

    #[test]
    fn mul_2_3() {
        let a = Fq::new(2.into());
        let b = Fq::new(3.into());

        let circuit = MulCircuit { a, b };

        let cs = ConstraintSystem::new_ref();
        circuit.generate_constraints(cs.clone()).unwrap();

        let is_satisfied = cs.is_satisfied().unwrap();
        if !is_satisfied {
            println!("{:?}", cs.which_is_unsatisfied().unwrap().unwrap());
        }
        arcworks_cs_to_pinocchio_r1cs(&cs);
        // con -
        // CS A: [[(Fp256(BigInteger256([10157024534604021774, 16668528035959406606, 5322190058819395602, 387181115924875961])), 1)]]

        //CS A: [[(Fp256(BigInteger256([9015221291577245683, 8239323489949974514, 1646089257421115374, 958099254763297437])), 1)]]
    }
}
