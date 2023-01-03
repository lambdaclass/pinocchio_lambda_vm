use crate::math::{self, group::Group};
use math::field_element::FieldElement as FE;

use super::prover::Proof;
use super::{prover::msm, setup::VerifyingKey};

pub fn verify(verifying_key: &VerifyingKey, proof: &Proof, c_input_output: &[FE]) -> bool {
    let b1 = check_divisibility(verifying_key, proof, c_input_output);
    let b2 = check_appropiate_spans(verifying_key, proof);
    let b3 = check_same_linear_combinations(verifying_key, proof);
    b1 && b2 && b3
}

pub fn check_divisibility(
    verifying_key: &VerifyingKey,
    proof: &Proof,
    input_output: &[FE],
) -> bool {
    let vk = verifying_key;

    let hiding_v = vk.gv_ks[0]
        .mul_by_group_element(msm(input_output, &vk.gv_ks[1..]))
        .mul_by_group_element(proof.g_vs);
    let hiding_w = vk.gw_ks[0]
        .mul_by_group_element(msm(input_output, &vk.gw_ks[1..]))
        .mul_by_group_element(proof.g_ws);
    let hiding_y = vk.gy_ks[0]
        .mul_by_group_element(msm(input_output, &vk.gy_ks[1..]))
        .mul_by_group_element(proof.g_ys);

    let lhs = hiding_v.pairing(hiding_w);
    let rhs_1 = verifying_key.gy_target_on_s.pairing(proof.g_hs);
    let rhs_2 = hiding_y.pairing(vk.g_1);

    lhs == rhs_1.mul_by_group_element(rhs_2)
}

pub fn check_appropiate_spans(verifying_key: &VerifyingKey, proof: &Proof) -> bool {
    let vk = verifying_key;

    proof.g_alpha_vs.pairing(vk.g_1) == proof.g_vs.pairing(vk.g_alpha_v)
        && proof.g_alpha_ws.pairing(vk.g_1) == proof.g_ws.pairing(vk.g_alpha_w)
        && proof.g_alpha_ys.pairing(vk.g_1) == proof.g_ys.pairing(vk.g_alpha_y)
}

pub fn check_same_linear_combinations(verifying_key: &VerifyingKey, proof: &Proof) -> bool {
    let vk = verifying_key;

    proof.g_beta_vwy.pairing(vk.g_gamma)
        == (proof
            .g_vs
            .mul_by_group_element(proof.g_ws)
            .mul_by_group_element(proof.g_ys))
        .pairing(vk.g_beta_gamma)
}
