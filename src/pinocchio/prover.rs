use crate::{
    circuits::qap::Qap,
    math::{self, group::Group},
};
use math::field_element::FieldElement as FE;

use super::setup::EvaluationKey;
pub type GroupType = FE;

/// Proof for Pinocchio
/// All but hs are the mid related elements
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Proof {
    g_vs: GroupType,
    g_ws: GroupType,
    g_ys: GroupType,
    g_alpha_vs: GroupType,
    g_alpha_ws: GroupType,
    g_alpha_ys: GroupType,
    g_beta_vwy: GroupType,
    //h may be the nil polynomial
    //and so g_hs will be empty
    g_hs: Option<GroupType>,
}

pub fn msm(c: &[FE], hidings: &[GroupType]) -> Option<GroupType> {
    c.iter()
        .zip(hidings.iter())
        .map(|(&c, &h)| h.mul_by_scalar(c))
        .reduce(|acc, x| acc + x)
}

/// Generates a proof from an evaluation_key,
/// the qap representation of the circuit
/// the first Cs of the inputs, the intermediate Cs results and the final Cs of the results
pub fn generate_proof(
    evaluation_key: &EvaluationKey,
    qap: &Qap,
    qap_c_coefficients: &[FE],
) -> Proof {
    let c_mid = &qap_c_coefficients
        [qap.number_of_inputs..(qap_c_coefficients.len() - qap.number_of_outputs)];

    let h_polynomial = qap.h_polynomial(qap_c_coefficients);

    // If C and the evaluation key are valid, this shouldn't fail
    // TO DO: Raise an error if some value can't be computated
    Proof {
        g_vs: msm(c_mid, &evaluation_key.gv_ks).unwrap(),
        g_ws: msm(c_mid, &evaluation_key.gw_ks).unwrap(),
        g_ys: msm(c_mid, &evaluation_key.gy_ks).unwrap(),
        g_alpha_vs: msm(c_mid, &evaluation_key.gv_alphaks).unwrap(),
        g_alpha_ws: msm(c_mid, &evaluation_key.gw_alphaks).unwrap(),
        g_alpha_ys: msm(c_mid, &evaluation_key.gy_alphaks).unwrap(),
        g_beta_vwy: msm(c_mid, &evaluation_key.g_beta).unwrap(),
        // evaluation key may have extra elements, since in this implementation
        // we used an upper bound of elements
        g_hs: msm(
            h_polynomial.coefficients(),
            &evaluation_key.g_s_i[..h_polynomial.coefficients().len()],
        ),
    }
}

#[cfg(test)]
mod tests {
    use crate::math::polynomial::Polynomial;

    use super::*;

    //MSM tests require the GroupType to be a FieldElement
    #[test]
    fn msm_11_is_1() {
        let c = [FE::one()];
        let hiding = [FE::one()];
        assert_eq!(msm(&c, &hiding).unwrap(), FE::one());
    }

    #[test]
    fn msm_23_is_6() {
        let c = [FE::new(3).unwrap()];
        let hiding = [FE::new(2).unwrap()];
        assert_eq!(msm(&c, &hiding).unwrap(), FE::new(6).unwrap());
    }

    #[test]
    fn msm_with_c_2_3_hiding_3_4_is_18() {
        let c = [FE::new(2).unwrap(), FE::new(3).unwrap()];
        let hiding = [FE::new(3).unwrap(), FE::new(4).unwrap()];
        assert_eq!(msm(&c, &hiding).unwrap(), FE::new(18).unwrap());
    }

    #[test]
    fn msm_with_empty_c_is_none() {
        let c = [];
        let hiding = [FE::new(3).unwrap(), FE::new(4).unwrap()];
        assert_eq!(msm(&c, &hiding), None);
    }

    #[test]
    fn msm_with_emtpy_hiding_is_none() {
        let c = [FE::zero()];
        let hiding = [];
        assert_eq!(msm(&c, &hiding), None);
    }

    #[test]
    fn msm_with_empty_arguments_is_none() {
        let c = [];
        let hiding = [];
        assert_eq!(msm(&c, &hiding), None);
    }

    // This test runs the proof algorithms with some easy inputs
    // to check operations are made correctly
    // Eval key and polynomials don't mean anything
    // It only works with FE as the GroupType
    #[test]
    fn proof_test() {
        // This eval key is the identity
        let easy_eval_key = EvaluationKey {
            gv_ks: vec![FE::one(), FE::one()],
            gw_ks: vec![FE::one(), FE::one()],
            gy_ks: vec![FE::one(), FE::one()],
            gv_alphaks: vec![FE::new(2).unwrap(), FE::new(2).unwrap()],
            gw_alphaks: vec![FE::new(2).unwrap(), FE::new(2).unwrap()],
            gy_alphaks: vec![FE::new(2).unwrap(), FE::new(2).unwrap()],
            g_s_i: vec![FE::one(), FE::one()],
            g_beta: vec![FE::one(), FE::one()],
        };

        // vwy are all equals, the 0 element is 0, the rest are ones
        let vw_polynomial = vec![
            // x0
            Polynomial::new(vec![FE::zero()]),
            //x1 = xinput
            Polynomial::new(vec![FE::one(), FE::one()]),
            //xmid
            Polynomial::new(vec![FE::one(), FE::one()]),
            Polynomial::new(vec![FE::one(), FE::one()]),
            //xoutput
            Polynomial::new(vec![FE::one(), FE::one()]),
        ];

        // y is different than vw, but only in values not in the middle
        // so results shouldn't change 
        let y_polynomial = vec![
            // x0
            Polynomial::new(vec![FE::new(32).unwrap()]),
            //x1 = xinput
            Polynomial::new(vec![FE::new(123).unwrap(), FE::new(123).unwrap()]),
            //xmid
            Polynomial::new(vec![FE::one(), FE::one()]),
            Polynomial::new(vec![FE::one(), FE::one()]),
            //xoutput
            Polynomial::new(vec![FE::new(321).unwrap(),FE::new(321).unwrap()]),
        ];

        // 1 input, and 1 output, so there are 2 values in the middle
        let easy_qap = Qap {
            v: vw_polynomial.clone(),
            w: vw_polynomial,
            y: y_polynomial,
            target: Polynomial::new(vec![FE::one(), FE::one()]),
            number_of_inputs: 1,
            number_of_outputs: 1,
        };

        let c_coefficients = vec![
            //c1
            FE::one(),
            //c_mids
            FE::new(2).unwrap(),
            FE::new(3).unwrap(),
            //c4 = c_output
            FE::one(),
        ];

        let proof = generate_proof(&easy_eval_key, &easy_qap, &c_coefficients);

        assert_eq!(proof.g_vs, FE::new(5).unwrap());
        assert_eq!(proof.g_ys, FE::new(5).unwrap());
        assert_eq!(proof.g_alpha_vs, FE::new(10).unwrap());
    }
}
