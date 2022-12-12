mod math;
mod pinocchio;
use math::field_element::FieldElement;
use pinocchio::proof::ToxicWaste;

fn main() {

    //should receive a circuit
    /*
    let proof = Prover::generate_proof(ToxicWaste::new());

    let verifier = Verifier::verify(proof);

    proof = Proof::new(ToxicWaste::new());

 */
    let toxic_waste = ToxicWaste::new();
}
