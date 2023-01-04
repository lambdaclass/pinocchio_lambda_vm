use crate::circuits::qap::Qap;
use crate::math::field_element::FieldElement as FE;
use crate::math::msm::msm;

use super::setup::EvaluationKey;
pub type CyclicGroupType = FE;

/// Proof for Pinocchio
/// All but hs are the mid related elements
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Proof {
    pub g_vs: CyclicGroupType,
    pub g_ws: CyclicGroupType,
    pub g_ys: CyclicGroupType,
    pub g_alpha_vs: CyclicGroupType,
    pub g_alpha_ws: CyclicGroupType,
    pub g_alpha_ys: CyclicGroupType,
    pub g_beta_vwy: CyclicGroupType,
    pub g_hs: CyclicGroupType,
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
        g_vs: msm(c_mid, &evaluation_key.gv_ks),
        g_ws: msm(c_mid, &evaluation_key.gw_ks),
        g_ys: msm(c_mid, &evaluation_key.gy_ks),
        g_alpha_vs: msm(c_mid, &evaluation_key.gv_alphaks),
        g_alpha_ws: msm(c_mid, &evaluation_key.gw_alphaks),
        g_alpha_ys: msm(c_mid, &evaluation_key.gy_alphaks),
        g_beta_vwy: msm(c_mid, &evaluation_key.g_beta),
        // evaluation key may have extra elements for g s i, since in this implementation it uses an upper bound of elements
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

    // This test runs the proof algorithms with some easy inputs
    // to check operations are made correctly
    // Eval key and polynomials don't mean anything
    // It only works with FE as the CyclicGroupType
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
            Polynomial::new(vec![FE::new(321).unwrap(), FE::new(321).unwrap()]),
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
