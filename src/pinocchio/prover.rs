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
    g_v_s: GroupType,
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
/// and the C of the C0, the inputs, the intermediate results and the final results
pub fn generate_proof(
    evaluation_key: &EvaluationKey,
    qap: &Qap,
    qap_c_coefficients: &[FE],
) -> Proof {
    //Inputs + solution should be an argument

    let c_mid = &qap_c_coefficients
        [qap.number_of_inputs + 1..(qap_c_coefficients.len() - qap.number_of_outputs)];

    let h_polynomial = qap.h_polynomial(qap_c_coefficients);

    // If C and the evaluation key are valid, this shouldn't fail
    // TO DO: Raise an error if some value can't be computated
    Proof {
        g_v_s: msm(c_mid, &evaluation_key.gv_ks).unwrap(),
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
