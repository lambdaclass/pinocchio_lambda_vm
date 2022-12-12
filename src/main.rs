mod math;
use math::field_element::FieldElement as FE;
use math::polynomial::Polynomial;

fn example_circuit() -> (FE, FE, Vec<Polynomial>, Vec<Polynomial>, Vec<Polynomial>) {
    let r5 = FE::new(0).unwrap();
    let r6 = FE::new(1).unwrap();

    let vs = vec![
        Polynomial::new(vec![r5, r6], vec![FE::zero(), FE::one()]),
        Polynomial::new(vec![r5, r6], vec![FE::zero(), FE::one()]),
        Polynomial::new(vec![r5, r6], vec![FE::one(), FE::zero()]),
        Polynomial::new(vec![r5, r6], vec![FE::zero(), FE::zero()]),
        Polynomial::new(vec![r5, r6], vec![FE::zero(), FE::zero()]),
        Polynomial::new(vec![r5, r6], vec![FE::zero(), FE::zero()]),
    ];

    let ws = vec![
        Polynomial::new(vec![r5, r6], vec![FE::zero(), FE::zero()]),
        Polynomial::new(vec![r5, r6], vec![FE::zero(), FE::zero()]),
        Polynomial::new(vec![r5, r6], vec![FE::zero(), FE::zero()]),
        Polynomial::new(vec![r5, r6], vec![FE::one(), FE::zero()]),
        Polynomial::new(vec![r5, r6], vec![FE::zero(), FE::one()]),
        Polynomial::new(vec![r5, r6], vec![FE::zero(), FE::zero()]),
    ];

    let ys = vec![
        Polynomial::new(vec![r5, r6], vec![FE::zero(), FE::zero()]),
        Polynomial::new(vec![r5, r6], vec![FE::zero(), FE::zero()]),
        Polynomial::new(vec![r5, r6], vec![FE::zero(), FE::zero()]),
        Polynomial::new(vec![r5, r6], vec![FE::zero(), FE::zero()]),
        Polynomial::new(vec![r5, r6], vec![FE::one(), FE::zero()]),
        Polynomial::new(vec![r5, r6], vec![FE::zero(), FE::one()]),
    ];

    (r5, r6, vs, ws, ys)
}

fn main() {
    // Example taken from the original Pinocchio paper at https://eprint.iacr.org/2013/279.pdf
    let (r5, r6, vs, ws, ys) = example_circuit();
    println!("    r5 r6");
    println!("---------");
    for (i, v_i) in vs.iter().enumerate() {
        println!("V_{:} {:}  {:}", i, v_i.evaluate(r5), v_i.evaluate(r6));
    }
    println!("---------");
    for (i, w_i) in ws.iter().enumerate() {
        println!("W_{:} {:}  {:}", i, w_i.evaluate(r5), w_i.evaluate(r6));
    }
    println!("---------");
    for (i, y_i) in ys.iter().enumerate() {
        println!("Y_{:} {:}  {:}", i, y_i.evaluate(r5), y_i.evaluate(r6));
    }
}
