use super::super::config::ORDER_R;
use super::prover::Proof;
use super::setup::VerifyingKey;
use crate::math::{self, cyclic_group::CyclicBilinearGroup};
use math::field_element::FieldElement;
use math::msm::msm;

type FE = FieldElement<ORDER_R>;

/// Pinocchio's verification algorithm.
pub fn verify<T: CyclicBilinearGroup>(
    verifying_key: &VerifyingKey<T>,
    proof: &Proof<T>,
    c_input_output: &[FE],
) -> bool {
    let b1 = check_divisibility(verifying_key, proof, c_input_output);
    let b2 = check_appropiate_spans(verifying_key, proof);
    let b3 = check_same_linear_combinations(verifying_key, proof);
    b1 && b2 && b3
}

pub fn check_divisibility<T: CyclicBilinearGroup>(
    verifying_key: &VerifyingKey<T>,
    proof: &Proof<T>,
    input_output: &[FE],
) -> bool {
    let vk = verifying_key;

    let hiding_v = vk.gv_ks[0]
        .operate_with(&msm(input_output, &vk.gv_ks[1..]))
        .operate_with(&proof.g_vs);
    let hiding_w = vk.gw_ks[0]
        .operate_with(&msm(input_output, &vk.gw_ks[1..]))
        .operate_with(&proof.g_ws);
    let hiding_y = vk.gy_ks[0]
        .operate_with(&msm(input_output, &vk.gy_ks[1..]))
        .operate_with(&proof.g_ys);

    let lhs = hiding_v.pairing(&hiding_w);
    let rhs_1 = verifying_key.gy_target_on_s.pairing(&proof.g_hs);
    let rhs_2 = hiding_y.pairing(&vk.g_1);

    lhs == rhs_1 * rhs_2
}

pub fn check_appropiate_spans<T: CyclicBilinearGroup>(
    verifying_key: &VerifyingKey<T>,
    proof: &Proof<T>,
) -> bool {
    let vk = verifying_key;

    let b1 = proof.g_alpha_vs.pairing(&vk.g_1) == proof.g_vs.pairing(&vk.g_alpha_v);
    let b2 = proof.g_alpha_ws.pairing(&vk.g_1) == proof.g_ws.pairing(&vk.g_alpha_w);
    let b3 = proof.g_alpha_ys.pairing(&vk.g_1) == proof.g_ys.pairing(&vk.g_alpha_y);
    b1 && b2 && b3
}

pub fn check_same_linear_combinations<T: CyclicBilinearGroup>(
    verifying_key: &VerifyingKey<T>,
    proof: &Proof<T>,
) -> bool {
    let vk = verifying_key;

    proof.g_beta_vwy.pairing(&vk.g_gamma)
        == (proof
            .g_vs
            .operate_with(&proof.g_ws)
            .operate_with(&proof.g_ys))
        .pairing(&vk.g_beta_gamma)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::elliptic_curve::EllipticCurveElement;

    // In this tests we do not hide elements. We work with the raw values instead.
    // These are easier to handle and computations can be done with pen and paper.
    // Check out the prover tests to see where most of these values come from.
    fn dummy_verifying_data_without_hidings() -> (VerifyingKey<FE>, Proof<FE>, Vec<FE>) {
        // Dummy verifying key assuming
        // (s, r_v, r_w, alpha_v, alpha_w, alpha_y, beta, gamma) = (1, 1, 1, 2, 2, 2, 3, 1)
        let verifying_key = VerifyingKey {
            g_1: FE::new(1),
            g_alpha_v: FE::new(2),
            g_alpha_w: FE::new(2),
            g_alpha_y: FE::new(2),
            g_gamma: FE::new(3),
            g_beta_gamma: FE::new(3),
            gy_target_on_s: FE::new(2),
            gv_ks: vec![FE::new(0), FE::new(6), FE::new(6)],
            gw_ks: vec![FE::new(0), FE::new(6), FE::new(6)],
            gy_ks: vec![FE::new(0), FE::new(6), FE::new(6)],
        };

        let input_output = vec![
            // One input value
            FE::new(3),
            // One output value
            FE::new(1),
        ];

        // This is a valid proof for the above input and output values.
        // See the prover tests to see where this comes from.
        let proof = Proof {
            g_vs: FE::new(6),
            g_ws: FE::new(6),
            g_ys: FE::new(6),
            g_hs: FE::new(120),
            g_alpha_vs: FE::new(12),
            g_alpha_ws: FE::new(12),
            g_alpha_ys: FE::new(12),
            g_beta_vwy: FE::new(18),
        };

        (verifying_key, proof, input_output)
    }

    #[test]
    fn test_divisibility_check_on_correct_proof_without_hidings() {
        let (verifying_key, proof, input_output) = dummy_verifying_data_without_hidings();
        assert!(check_divisibility(&verifying_key, &proof, &input_output));
    }

    #[test]
    fn test_appropiate_spans_on_correct_proof_without_hidings() {
        let (verifying_key, proof, _) = dummy_verifying_data_without_hidings();
        assert!(check_appropiate_spans(&verifying_key, &proof));
    }

    #[test]
    fn test_same_linear_combinations_on_correct_proof_without_hidings() {
        let (verifying_key, proof, _) = dummy_verifying_data_without_hidings();
        assert!(check_same_linear_combinations(&verifying_key, &proof));
    }

    #[test]
    fn test_divisibility_check_correct_on_incorrect_input_output_without_hidings() {
        let (verifying_key, proof, _) = dummy_verifying_data_without_hidings();
        let input_output = vec![FE::new(0), FE::new(0)];
        assert!(!check_divisibility(&verifying_key, &proof, &input_output));
    }

    #[test]
    fn test_appropiate_spans_correct_on_incorrect_proof_without_hidings() {
        let (verifying_key, mut proof, _) = dummy_verifying_data_without_hidings();
        proof.g_alpha_vs = FE::new(0);
        assert!(!check_appropiate_spans(&verifying_key, &proof));
    }

    #[test]
    fn test_same_linear_combinations_on_incorrect_proof_without_hidings() {
        let (verifying_key, mut proof, _) = dummy_verifying_data_without_hidings();
        proof.g_beta_vwy = FE::new(0);
        assert!(!check_same_linear_combinations(&verifying_key, &proof));
    }

    // The following is the same as before but with hidings on elliptic curves.
    fn dummy_verifying_data_with_elliptic_curve_hidings() -> (
        VerifyingKey<EllipticCurveElement>,
        Proof<EllipticCurveElement>,
        Vec<FE>,
    ) {
        let g = EllipticCurveElement::generator();
        // Dummy verifying key assuming
        // (s, r_v, r_w, alpha_v, alpha_w, alpha_y, beta, gamma) = (1, 1, 1, 2, 2, 2, 3, 1)
        let verifying_key = VerifyingKey {
            g_1: g.operate_with_self(1),
            g_alpha_v: g.operate_with_self(2),
            g_alpha_w: g.operate_with_self(2),
            g_alpha_y: g.operate_with_self(2),
            g_gamma: g.operate_with_self(3),
            g_beta_gamma: g.operate_with_self(3),
            gy_target_on_s: g.operate_with_self(2),
            gv_ks: vec![
                g.operate_with_self(0),
                g.operate_with_self(6),
                g.operate_with_self(6),
            ],
            gw_ks: vec![
                g.operate_with_self(0),
                g.operate_with_self(6),
                g.operate_with_self(6),
            ],
            gy_ks: vec![
                g.operate_with_self(0),
                g.operate_with_self(6),
                g.operate_with_self(6),
            ],
        };

        let input_output = vec![
            // One input value
            FE::new(3),
            // One output value
            FE::new(1),
        ];

        // This is a valid proof for the above input and output values.
        // See the prover tests to see where this comes from.
        let proof = Proof {
            g_vs: g.operate_with_self(6),
            g_ws: g.operate_with_self(6),
            g_ys: g.operate_with_self(6),
            g_hs: g.operate_with_self(120),
            g_alpha_vs: g.operate_with_self(12),
            g_alpha_ws: g.operate_with_self(12),
            g_alpha_ys: g.operate_with_self(12),
            g_beta_vwy: g.operate_with_self(18),
        };

        (verifying_key, proof, input_output)
    }

    #[test]
    fn test_divisibility_check_correct_on_correct_proof_with_elliptic_curve_hidings() {
        let (verifying_key, proof, input_output) =
            dummy_verifying_data_with_elliptic_curve_hidings();
        assert!(check_divisibility(&verifying_key, &proof, &input_output));
    }

    #[test]
    fn test_appropiate_spans_correct_on_correct_proof_with_elliptic_curve_hidings() {
        let (verifying_key, proof, _) = dummy_verifying_data_with_elliptic_curve_hidings();
        assert!(check_appropiate_spans(&verifying_key, &proof));
    }

    #[test]
    fn test_same_linear_combinations_on_correct_proof_with_elliptic_curve_hidings() {
        let (verifying_key, proof, _) = dummy_verifying_data_with_elliptic_curve_hidings();
        assert!(check_same_linear_combinations(&verifying_key, &proof));
    }

    #[test]
    fn test_divisibility_check_correct_on_incorrect_input_output_with_elliptic_curve_hidings() {
        let (verifying_key, proof, _) = dummy_verifying_data_with_elliptic_curve_hidings();
        let input_output = vec![FE::new(0), FE::new(1)];
        assert!(!check_divisibility(&verifying_key, &proof, &input_output));
    }

    #[test]
    fn test_appropiate_spans_correct_on_incorrect_proof_with_elliptic_curve_hidings() {
        let (verifying_key, mut proof, _) = dummy_verifying_data_with_elliptic_curve_hidings();
        proof.g_alpha_vs = EllipticCurveElement::neutral_element();
        assert!(!check_appropiate_spans(&verifying_key, &proof));
    }

    #[test]
    fn test_same_linear_combinations_on_incorrect_proof_with_elliptic_curve_hidings() {
        let (verifying_key, mut proof, _) = dummy_verifying_data_with_elliptic_curve_hidings();
        proof.g_beta_vwy = EllipticCurveElement::neutral_element();
        assert!(!check_same_linear_combinations(&verifying_key, &proof));
    }
}
