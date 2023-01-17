use super::super::config::ORDER_R;
use crate::circuits::qap::QuadraticArithmeticProgram as QAP;
use crate::math::cyclic_group::CyclicBilinearGroup;
use crate::math::field_element::FieldElement;
use crate::math::msm::msm;

use super::setup::EvaluationKey;

type FE = FieldElement<ORDER_R>;

/// Pinocchio's proof
/// Using the notation of Pinocchio's paper, these are
/// the hidings of v_{mid}(s), w_{mid}(s), y_{mid}(s), h(s)
/// and the "redundant" hidings for the consistency checks
/// of the verifier.
/// The polynomials v, w and y are the polynomials from the QAP
/// The polynomial h is equal to (vw - y) / t, where t is the
/// target polynomial of the QAP.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Proof<T: CyclicBilinearGroup> {
    pub g_vs: T,
    pub g_ws: T,
    pub g_ys: T,
    pub g_hs: T,
    pub g_alpha_vs: T,
    pub g_alpha_ws: T,
    pub g_alpha_ys: T,
    pub g_beta_vwy: T,
}

/// Generates a proof.
/// Takes as input the evaluation_key,
/// the QAP representation of the circuit
/// and the values of the circuit wires corresponding
/// to the particular execution instance. These values are
/// the ones denoted `c_i` in the paper. They include all
/// inputs and outputs.
pub fn generate_proof<T: CyclicBilinearGroup>(
    evaluation_key: &EvaluationKey<T>,
    qap: &QAP,
    qap_c_coefficients: &[FE],
) -> Proof<T> {
    let c_mid = &qap_c_coefficients
        [qap.number_of_inputs..(qap_c_coefficients.len() - qap.number_of_outputs)];

    let h_polynomial = qap.h_polynomial(qap_c_coefficients);

    Proof {
        g_vs: msm(c_mid, &evaluation_key.gv_ks),
        g_ws: msm(c_mid, &evaluation_key.gw_ks),
        g_ys: msm(c_mid, &evaluation_key.gy_ks),
        g_alpha_vs: msm(c_mid, &evaluation_key.gv_alphaks),
        g_alpha_ws: msm(c_mid, &evaluation_key.gw_alphaks),
        g_alpha_ys: msm(c_mid, &evaluation_key.gy_alphaks),
        g_beta_vwy: msm(c_mid, &evaluation_key.g_beta),
        g_hs: msm(
            h_polynomial.coefficients(),
            &evaluation_key.g_s_i[..h_polynomial.coefficients().len()],
        ),
    }
}

#[cfg(test)]
mod tests {
    use crate::math::elliptic_curve::EllipticCurveElement;
    use crate::math::{cyclic_group::CyclicBilinearGroup, polynomial::Polynomial};

    use super::*;

    // This test runs the proof algorithms with some easy inputs
    // to check operations are correct
    // The evaluation key and polynomials have dummy values and do
    // not correspond to any circuit.
    #[test]
    fn test_proof_is_generated_correctly_when_using_raw_values_instead_of_hidings() {
        // In this test we do not hide elements. We work with the raw values instead.
        // These are easier to handle and computations can be done with pen and paper.

        // We choose a dummy `target_qap` to be t = X + 1
        let target_qap = Polynomial::new(vec![FE::new(1), FE::new(1)]);

        // We choose vs = [t, t, t, t, t]. The same for ws, and ys.
        let vs = vec![target_qap.clone(); 5];
        let ws = vec![target_qap.clone(); 5];
        let ys = vec![target_qap.clone(); 5];

        // There is 1 input and 1 output. So there are 2 middle values.
        let easy_qap = QAP {
            vs,
            ws,
            ys,
            target: target_qap,
            number_of_inputs: 1,
            number_of_outputs: 1,
        };
        // Dummy evaluation key assuming
        // (s, r_v, r_w, alpha_v, alpha_w, alpha_y, beta, gamma) = (1, 1, 1, 2, 2, 2, 3, 1)
        let evaluation_key = EvaluationKey {
            gv_ks: vec![FE::new(2), FE::new(2)],
            gw_ks: vec![FE::new(2), FE::new(2)],
            gy_ks: vec![FE::new(2), FE::new(2)],
            gv_alphaks: vec![FE::new(4), FE::new(4)],
            gw_alphaks: vec![FE::new(4), FE::new(4)],
            gy_alphaks: vec![FE::new(4), FE::new(4)],
            g_s_i: vec![FE::new(1), FE::new(1)],
            g_beta: vec![FE::new(18), FE::new(18)],
        };

        let c_coefficients = vec![
            // c1
            FE::new(3),
            // c_mids:
            // c_2
            FE::new(2),
            // c_3
            FE::new(1),
            // output:
            // c4
            FE::new(1),
        ];

        let proof = generate_proof(&evaluation_key, &easy_qap, &c_coefficients);

        // For this choice of dummy polynomials and circuit wire values, we obtain
        // v = t + c_1 * t + c_2 * t + c_3 * t + c_4 * t, and v_{mid} = c_2 * t + c_3 * t.
        // And similarly for v and w.
        // Since in our evaluation key we are assuming that the component `s`
        // of the toxic waste is 1, we have v_{mid}(s) = 2 * t(1) + 1 * t(1) = 6.
        assert_eq!(proof.g_vs, FE::new(6));
        assert_eq!(proof.g_ws, FE::new(6));
        assert_eq!(proof.g_ys, FE::new(6));
        assert_eq!(proof.g_alpha_vs, FE::new(2 * 6));
        assert_eq!(proof.g_alpha_ws, FE::new(2 * 6));
        assert_eq!(proof.g_alpha_ws, FE::new(2 * 6));

        // On the other hand p = vw - y = r^2 * t^2 - r * t,
        // where r = 1 + c_1 + c_2 + c_3 + c_4 = 8.
        // Therefore h = p / t = r^2 * t - r = 64 * t - 8.
        // we have h(s) = h(1) = 64 * t(1) - 8 = 64 * 2 - 8 = 120.
        assert_eq!(proof.g_hs, FE::new(120));

        assert_eq!(proof.g_beta_vwy, FE::new(3 * 6 + 3 * 6 + 3 * 6));
    }

    #[test]
    fn test_proof_is_generated_correctly_when_hiding_in_elliptic_curves() {
        // Same test as before, but this time hiding points in elliptic curves
        let target_qap = Polynomial::new(vec![FE::new(1), FE::new(1)]);
        let vs = vec![target_qap.clone(); 5];
        let ws = vec![target_qap.clone(); 5];
        let ys = vec![target_qap.clone(); 5];

        let easy_qap = QAP {
            vs,
            ws,
            ys,
            target: target_qap,
            number_of_inputs: 1,
            number_of_outputs: 1,
        };

        let g = EllipticCurveElement::generator();

        let evaluation_key = EvaluationKey {
            gv_ks: vec![g.operate_with_self(2), g.operate_with_self(2)],
            gw_ks: vec![g.operate_with_self(2), g.operate_with_self(2)],
            gy_ks: vec![g.operate_with_self(2), g.operate_with_self(2)],
            gv_alphaks: vec![g.operate_with_self(4), g.operate_with_self(4)],
            gw_alphaks: vec![g.operate_with_self(4), g.operate_with_self(4)],
            gy_alphaks: vec![g.operate_with_self(4), g.operate_with_self(4)],
            g_s_i: vec![g.clone(), g.clone()],
            g_beta: vec![g.operate_with_self(18), g.operate_with_self(18)],
        };

        let c_coefficients = vec![FE::new(3), FE::new(2), FE::new(1), FE::new(1)];

        let proof = generate_proof(&evaluation_key, &easy_qap, &c_coefficients);

        assert_eq!(proof.g_vs, g.operate_with_self(6));
        assert_eq!(proof.g_ws, g.operate_with_self(6));
        assert_eq!(proof.g_ys, g.operate_with_self(6));
        assert_eq!(proof.g_alpha_vs, g.operate_with_self(2 * 6));
        assert_eq!(proof.g_alpha_ws, g.operate_with_self(2 * 6));
        assert_eq!(proof.g_alpha_ws, g.operate_with_self(2 * 6));

        assert_eq!(proof.g_hs, g.operate_with_self(120));

        assert_eq!(proof.g_beta_vwy, g.operate_with_self(3 * 6 + 3 * 6 + 3 * 6));
    }
}
