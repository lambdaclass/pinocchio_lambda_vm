use crate::circuits::qap::Qap;
use crate::math::field_element::FieldElement;
use crate::math::msm::msm;
use crate::math::elliptic_curve::EllipticCurveElement;

use super::setup::EvaluationKey;

const ORDER: u128 = 5;
type FE = FieldElement<ORDER>;
pub type CyclicGroupType = EllipticCurveElement;

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
    use crate::math::{polynomial::Polynomial, cyclic_group::CyclicGroup};

    use super::*;

    // This test runs the proof algorithms with some easy inputs
    // to check operations are made correctly
    // Eval key and polynomials don't mean anything
    // It only works with FE as the CyclicGroupType
    #[test]
    fn proof_test() {
        let g = CyclicGroupType::generator();

        // This eval key is the identity
        let easy_eval_key = EvaluationKey {
            gv_ks: vec![g.clone(), g.clone()],
            gw_ks: vec![g.clone(), g.clone()],
            gy_ks: vec![g.clone(), g.clone()],
            gv_alphaks: vec![g.operate_with_self(2), g.operate_with_self(2)],
            gw_alphaks: vec![g.operate_with_self(2), g.operate_with_self(2)],
            gy_alphaks: vec![g.operate_with_self(2), g.operate_with_self(2)],
            g_s_i: vec![g.clone(), g.clone()],
            g_beta: vec![g.clone(), g.clone()],
        };

        // vwy are all equals, the 0 element is 0, the rest are ones
        let vw_polynomial = vec![
            // x0
            Polynomial::new(vec![FE::new(0)]),
            //x1 = xinput
            Polynomial::new(vec![FE::new(1), FE::new(1)]),
            //xmid
            Polynomial::new(vec![FE::new(1), FE::new(1)]),
            Polynomial::new(vec![FE::new(1), FE::new(1)]),
            //xoutput
            Polynomial::new(vec![FE::new(1), FE::new(1)]),
        ];

        // y is different than vw, but only in values not in the middle
        // so results shouldn't change
        let y_polynomial = vec![
            // x0
            Polynomial::new(vec![FE::new(32)]),
            //x1 = xinput
            Polynomial::new(vec![FE::new(123), FE::new(123)]),
            //xmid
            Polynomial::new(vec![FE::new(1), FE::new(1)]),
            Polynomial::new(vec![FE::new(1), FE::new(1)]),
            //xoutput
            Polynomial::new(vec![FE::new(321), FE::new(321)]),
        ];

        // 1 input, and 1 output, so there are 2 values in the middle
        let easy_qap = Qap {
            v: vw_polynomial.clone(),
            w: vw_polynomial,
            y: y_polynomial,
            target: Polynomial::new(vec![FE::new(1), FE::new(1)]),
            number_of_inputs: 1,
            number_of_outputs: 1,
        };

        let c_coefficients = vec![
            //c1
            FE::new(1),
            //c_mids
            FE::new(2),
            FE::new(3),
            //c4 = c_output
            FE::new(1),
        ];

        let proof = generate_proof(&easy_eval_key, &easy_qap, &c_coefficients);

        assert_eq!(proof.g_vs, g.operate_with_self(5));
        assert_eq!(proof.g_ys, g.operate_with_self(5));
        assert_eq!(proof.g_alpha_vs, g.operate_with_self(5));
    }
}
