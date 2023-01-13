use ark_bn254::{Fq, FqParameters};
use ark_ff::{BigInteger256, Fp256, FpParameters, PrimeField};
use ark_relations::r1cs::ConstraintSystemRef;
use num_bigint::BigUint;
use pinocchio_vm::circuits::r1cs::R1CS;

pub fn arcworks_cs_to_pinocchio_r1cs(cs: &ConstraintSystemRef<Fp256<FqParameters>>) {
    cs.inline_all_lcs();

    let a: Vec<(String, usize)> = cs.to_matrices().unwrap().a[0]
        .iter()
        .map(|(x, y)| (x.to_string(), *y))
        .collect();
    let b: Vec<(String, usize)> = cs.to_matrices().unwrap().b[0]
        .iter()
        .map(|(x, y)| (x.to_string(), *y))
        .collect();
    let c: Vec<(String, usize)> = cs.to_matrices().unwrap().c[0]
        .iter()
        .map(|(x, y)| (x.to_string(), *y))
        .collect();

    println!("Mod: {:?}", FqParameters::MODULUS.to_string());
    println!("R: {:?}", FqParameters::R.to_string());
    println!("Inv: {:?}", FqParameters::INV.to_string());


    println!("CS A: {:?}", a);
    println!("CS B: {:?}", b);
    println!("CS B: {:?}", c);

    let r1cs_matrixes = cs.to_matrices().unwrap();

    let a = sparse_arcworks_matrix_fp_to_u128(&r1cs_matrixes.a);
    let b = sparse_arcworks_matrix_fp_to_u128(&r1cs_matrixes.b);
    let c = sparse_arcworks_matrix_fp_to_u128(&r1cs_matrixes.c);

    println!("Converted CS A: {:?}", a);
    println!("Converted CS B: {:?}", b);
    println!("Converted CS B: {:?}", c);
}

fn sparse_arcworks_matrix_fp_to_u128 (m: &Vec<Vec<(Fp256<FqParameters>,usize)>>) -> Vec<Vec<(u128, usize)>> 
    {
    m.iter().map( |x| {
        x.iter()
        .map( |(x, y)| {
            (biguint_to_u128(x.into_repr().into()),*y)
        })
        .collect()
    }).collect()
}


/// Converts biguint to u128
/// If the biguint is bigger than u128 it takes the first 128 bytes
fn biguint_to_u128(big: BigUint) -> u128 {

    match big.to_u64_digits().len() {
        0 => 0,
        1 => big.to_u64_digits()[0].into(),
        _ => 
            big.to_u64_digits()[0] as u128 
            & ((big.to_u64_digits()[1] as u128) << 64)
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
    fn mul_with_add() {
        let a = Fq::new(2.into());
        let b = Fq::new(3.into());

        let circuit = MulCircuit { a, b };

        let big: BigInteger256 = 7.into();
        println!("Big: {:?}", big);

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
