use crate::{circuits::qap::Qap, math};
use math::field_element::FieldElement as FE;

use super::setup::EvaluationKey;
pub type GroupType = FE;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Proof {
    g_vmid_s: GroupType,
    g_ws: GroupType,
    g_ys: GroupType,
    g_hs: GroupType,
    g_alpha_vmids: GroupType,
    g_alpha_ws: GroupType,
    g_alpha_ys: GroupType,
    g_alpha_hs: GroupType,
    g_beta_vwy: GroupType,
}

#[allow(dead_code)]
pub fn generate_proof(_proving_key: EvaluationKey, _qap: Qap) -> Proof {
    Proof {
        g_vmid_s: GroupType::zero(),
        g_ws: GroupType::zero(),
        g_ys: GroupType::zero(),
        g_hs: GroupType::zero(),
        g_alpha_vmids: GroupType::zero(),
        g_alpha_ws: GroupType::zero(),
        g_alpha_ys: GroupType::zero(),
        g_alpha_hs: GroupType::zero(),
        g_beta_vwy: GroupType::zero(),
    }
}
