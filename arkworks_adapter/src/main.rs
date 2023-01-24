use arkworks_adapter::{self, fq5, pinocchio_r1cs_from_arkworks_cs, pinocchio_io_and_witness_from_arkworks_cs};
use pinocchio_lambda_vm::{pinocchio::{prover, verifier, setup::{EvaluationKey, VerificationKey, setup, ToxicWaste}}, math::elliptic_curve::EllipticCurveElement, circuits::qap::QuadraticArithmeticProgram};
use poc_encryption_proof::{self, encrypt_cs};
use aes::{
    cipher::{BlockEncrypt, KeyInit, generic_array::GenericArray},
    Aes128,
};
fn main() {
    let message = [1_u8; 16];
    let secret_key = [0_u8; 16];

    let primitive_secret_key = Aes128::new(GenericArray::from_slice(&secret_key));

    let primitive_ciphertext = primitive_encrypt(&message, &primitive_secret_key);

    let cs = encrypt_cs::<fq5::Fq>(&message,&secret_key,&primitive_ciphertext).unwrap();

    println!("CS");
    println!("Satisfied: {:?}", cs.is_satisfied());
    println!("Constraints: {}", cs.num_constraints());
    println!("Instance variables: {}", cs.num_instance_variables());
    println!("Witness variables: {}", cs.num_witness_variables());

    let toxic_waste = ToxicWaste::sample();

    let r1cs = pinocchio_r1cs_from_arkworks_cs(&cs);
    let (io, witness) = pinocchio_io_and_witness_from_arkworks_cs(&cs);

    println!("Converted to Pinocchio format");
    println!("IO: {:?}", io);


    let qap: QuadraticArithmeticProgram = r1cs.into();

    println!("Generated QAP");
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

    println!("Verification passed?: {}", accepted);

    assert!(accepted);
}

pub fn primitive_encrypt(message: &[u8; 16], primitive_secret_key: &Aes128) -> Vec<u8> {
    let mut encrypted_message = Vec::new();
    let mut block = GenericArray::clone_from_slice(message);
    primitive_secret_key.encrypt_block(&mut block);
    encrypted_message.extend_from_slice(block.as_slice());
    encrypted_message
}
