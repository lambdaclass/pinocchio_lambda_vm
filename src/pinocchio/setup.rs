use super::super::config::ORDER_R;
use crate::circuits::qap::QuadraticArithmeticProgram as QAP;
use crate::math;
use math::cyclic_group::CyclicBilinearGroup;
use math::field_element::FieldElement;

pub type FE = FieldElement<ORDER_R>;

/// Evaluation key for Pinocchio
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EvaluationKey<T: CyclicBilinearGroup> {
    pub gv_ks: Vec<T>,
    pub gw_ks: Vec<T>,
    pub gy_ks: Vec<T>,
    pub gv_alphaks: Vec<T>,
    pub gw_alphaks: Vec<T>,
    pub gy_alphaks: Vec<T>,
    pub g_s_i: Vec<T>,
    pub g_beta: Vec<T>,
}
/// Verifying key for Pinocchio
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VerificationKey<T: CyclicBilinearGroup> {
    pub g_1: T,
    pub g_alpha_v: T,
    pub g_alpha_w: T,
    pub g_alpha_y: T,
    pub g_gamma: T,
    pub g_beta_gamma: T,
    pub gy_target_on_s: T,
    pub gv_ks: Vec<T>,
    pub gw_ks: Vec<T>,
    pub gy_ks: Vec<T>,
}
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct ToxicWaste {
    s: FE,
    alpha_v: FE,
    alpha_w: FE,
    alpha_y: FE,
    beta: FE,
    rv: FE,
    rw: FE,
    gamma: FE,
}

impl ToxicWaste {
    pub fn ry(self) -> FE {
        self.rv * self.rw
    }

    #[allow(clippy::too_many_arguments)]
    pub fn new(
        s: FE,
        alpha_v: FE,
        alpha_w: FE,
        alpha_y: FE,
        beta: FE,
        rv: FE,
        rw: FE,
        gamma: FE,
    ) -> Self {
        Self {
            s,
            alpha_v,
            alpha_w,
            alpha_y,
            beta,
            rv,
            rw,
            gamma,
        }
    }

    pub fn sample() -> Self {
        Self {
            s: FE::random(),
            alpha_v: FE::random(),
            alpha_w: FE::random(),
            alpha_y: FE::random(),
            beta: FE::random(),
            rv: FE::random(),
            rw: FE::random(),
            gamma: FE::random(),
        }
    }
}

fn generate_verification_key<T: CyclicBilinearGroup>(
    qap: &QAP,
    toxic_waste: &ToxicWaste,
    generator: &T,
) -> VerificationKey<T> {
    let s = toxic_waste.s;
    let alpha_v = toxic_waste.alpha_v;
    let alpha_w = toxic_waste.alpha_w;
    let alpha_y = toxic_waste.alpha_y;
    let beta = toxic_waste.beta;
    let rv = toxic_waste.rv;
    let rw = toxic_waste.rw;
    let gamma = toxic_waste.gamma;
    let ry = toxic_waste.ry();

    let g = generator;

    let vector_capacity = qap.number_of_inputs + qap.number_of_inputs + 1;
    let mut gv_ks_io: Vec<T> = Vec::with_capacity(vector_capacity);
    let mut gw_ks_io: Vec<T> = Vec::with_capacity(vector_capacity);
    let mut gy_ks_io: Vec<T> = Vec::with_capacity(vector_capacity);

    gv_ks_io.push(g.operate_with_self((rv * qap.v0().evaluate(s)).representative()));
    gw_ks_io.push(g.operate_with_self((rw * qap.w0().evaluate(s)).representative()));
    gy_ks_io.push(g.operate_with_self((ry * qap.y0().evaluate(s)).representative()));

    for k in 0..qap.v_input().len() {
        gv_ks_io.push(g.operate_with_self((rv * qap.v_input()[k].evaluate(s)).representative()));
        gw_ks_io.push(g.operate_with_self((rw * qap.w_input()[k].evaluate(s)).representative()));
        gy_ks_io.push(g.operate_with_self((ry * qap.y_input()[k].evaluate(s)).representative()));
    }

    for k in 0..qap.v_output().len() {
        gv_ks_io.push(g.operate_with_self((rv * qap.v_output()[k].evaluate(s)).representative()));
        gw_ks_io.push(g.operate_with_self((rw * qap.w_output()[k].evaluate(s)).representative()));
        gy_ks_io.push(g.operate_with_self((ry * qap.y_output()[k].evaluate(s)).representative()));
    }

    VerificationKey {
        g_1: g.clone(),
        g_alpha_v: g.operate_with_self(alpha_v.representative()),
        g_alpha_w: g.operate_with_self(alpha_w.representative()),
        g_alpha_y: g.operate_with_self(alpha_y.representative()),
        g_gamma: g.operate_with_self(gamma.representative()),
        g_beta_gamma: g.operate_with_self((beta * gamma).representative()),
        gy_target_on_s: g.operate_with_self((ry * qap.target.evaluate(s)).representative()),
        gv_ks: gv_ks_io,
        gw_ks: gw_ks_io,
        gy_ks: gy_ks_io,
    }
}

fn generate_evaluation_key<T: CyclicBilinearGroup>(
    qap: &QAP,
    toxic_waste: &ToxicWaste,
    generator: &T,
) -> EvaluationKey<T> {
    let (vs_mid, ws_mid, ys_mid) = (qap.v_mid(), qap.w_mid(), qap.y_mid());

    let s = toxic_waste.s;
    let alpha_v = toxic_waste.alpha_v;
    let alpha_w = toxic_waste.alpha_w;
    let alpha_y = toxic_waste.alpha_y;
    let beta = toxic_waste.beta;
    let rv = toxic_waste.rv;
    let rw = toxic_waste.rw;
    let ry = toxic_waste.ry();

    let g = generator;

    let degree = qap.target.degree();

    let mut gv_ks_mid: Vec<T> = Vec::with_capacity(vs_mid.len());
    let mut gw_ks_mid: Vec<T> = Vec::with_capacity(vs_mid.len());
    let mut gy_ks_mid: Vec<T> = Vec::with_capacity(vs_mid.len());
    let mut gv_alphaks_mid: Vec<T> = Vec::with_capacity(vs_mid.len());
    let mut gw_alphaks_mid: Vec<T> = Vec::with_capacity(vs_mid.len());
    let mut gy_alphaks_mid: Vec<T> = Vec::with_capacity(vs_mid.len());
    let mut g_beta_mid: Vec<T> = Vec::with_capacity(vs_mid.len());
    // g_s_i is the only paramater to depend on the degree of the qap
    // This is an upper bound, it could be smaller
    let mut g_s_i: Vec<T> = Vec::with_capacity(degree);

    // Set evaluation keys for each of their respective k mid element
    for k in 0..vs_mid.len() {
        gv_ks_mid.push(g.operate_with_self((rv * vs_mid[k].evaluate(s)).representative()));
        gw_ks_mid.push(g.operate_with_self((rw * ws_mid[k].evaluate(s)).representative()));
        gy_ks_mid.push(g.operate_with_self((ry * ys_mid[k].evaluate(s)).representative()));
        gv_alphaks_mid
            .push(g.operate_with_self((ry * alpha_v * vs_mid[k].evaluate(s)).representative()));
        gw_alphaks_mid
            .push(g.operate_with_self((rw * alpha_w * ws_mid[k].evaluate(s)).representative()));
        gy_alphaks_mid
            .push(g.operate_with_self((ry * alpha_y * ys_mid[k].evaluate(s)).representative()));
        g_beta_mid.push(
            g.operate_with_self(
                (rv * beta * vs_mid[k].evaluate(s)
                    + rw * beta * ws_mid[k].evaluate(s)
                    + ry * beta * ys_mid[k].evaluate(s))
                .representative(),
            ),
        )
    }

    for i in 0..qap.target.degree() {
        // This unwrap would only fail in an OS
        // with 256 bits pointer, which doesn't exist
        g_s_i.push(g.operate_with_self(s.pow(i.try_into().unwrap()).representative()));
    }

    EvaluationKey {
        gv_ks: gv_ks_mid,
        gw_ks: gw_ks_mid,
        gy_ks: gy_ks_mid,
        gv_alphaks: gv_alphaks_mid,
        gw_alphaks: gw_alphaks_mid,
        gy_alphaks: gy_alphaks_mid,
        g_s_i,
        g_beta: g_beta_mid,
    }
}

pub fn setup<T: CyclicBilinearGroup>(
    qap: &QAP,
    toxic_waste: &ToxicWaste,
) -> (EvaluationKey<T>, VerificationKey<T>) {
    let generator = T::generator();
    (
        generate_evaluation_key(qap, toxic_waste, &generator),
        generate_verification_key(qap, toxic_waste, &generator),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::{setup, ToxicWaste};
    use crate::circuits::test_utils::new_test_qap;

    fn identity_toxic_waste() -> ToxicWaste {
        ToxicWaste {
            s: FE::new(1),
            alpha_v: FE::new(1),
            alpha_w: FE::new(1),
            alpha_y: FE::new(1),
            beta: FE::new(1),
            rv: FE::new(1),
            rw: FE::new(1),
            gamma: FE::new(1),
        }
    }

    #[test]
    fn evaluation_keys_size_for_test_circuit_is_1_for_each_key() {
        let (eval_key, _): (EvaluationKey<FE>, VerificationKey<FE>) =
            setup(&new_test_qap(), &identity_toxic_waste());
        assert_eq!(eval_key.gv_ks.len(), 1);
        assert_eq!(eval_key.gw_ks.len(), 1);
        assert_eq!(eval_key.gy_ks.len(), 1);
        assert_eq!(eval_key.gv_alphaks.len(), 1);
        assert_eq!(eval_key.gw_alphaks.len(), 1);
        assert_eq!(eval_key.gy_alphaks.len(), 1);
        assert_eq!(eval_key.g_beta.len(), 1);
    }

    // This test makes a manual computation of each field and compares it
    // with the given values of the eval key
    #[test]
    fn eval_key_returns_appropiate_values() {
        let r5 = FE::new(2);
        let tw = ToxicWaste {
            s: FE::new(2),
            alpha_v: FE::new(2),
            alpha_w: FE::new(2),
            alpha_y: FE::new(2),
            beta: FE::new(2),
            rv: FE::new(2),
            rw: FE::new(2),
            gamma: FE::new(1),
        };

        let g = FE::generator();
        let test_circuit = new_test_qap();

        let (eval_key, _): (EvaluationKey<FE>, VerificationKey<FE>) = setup(&test_circuit, &tw);

        // These keys should be the same evaluation * rv, which is two
        assert_eq!(
            eval_key.gv_ks[0],
            g.operate_with_self(
                (test_circuit.v_mid()[0].evaluate(r5) * FE::new(2)).representative()
            )
        );
        assert_eq!(
            eval_key.gw_ks[0],
            g.operate_with_self(
                (test_circuit.w_mid()[0].evaluate(r5) * FE::new(2)).representative()
            )
        );
        // These keys should be the same evaluation * ys, which is two
        // Since the whole thing is 0
        assert_eq!(
            eval_key.gy_ks[0],
            g.operate_with_self(
                (test_circuit.y_mid()[0].evaluate(r5) * FE::new(4)).representative()
            )
        );

        // alpha * rv and alpha * rw is 4
        assert_eq!(
            eval_key.gv_alphaks[0],
            g.operate_with_self(
                (test_circuit.v_mid()[0].evaluate(r5) * FE::new(4)).representative()
            )
        );
        assert_eq!(
            eval_key.gv_alphaks[0],
            g.operate_with_self(
                (test_circuit.v_mid()[0].evaluate(r5) * FE::new(4)).representative()
            )
        );
        // alpha * ry and alpha * rw is 8
        assert_eq!(
            eval_key.gv_alphaks[0],
            g.operate_with_self(
                (test_circuit.v_mid()[0].evaluate(r5) * FE::new(8)).representative()
            )
        );

        assert_eq!(
            eval_key.g_beta[0],
            // beta * rv is 4
            g.operate_with_self(
                (test_circuit.v_mid()[0].evaluate(r5) * FE::new(4) +
            test_circuit.w_mid()[0].evaluate(r5) * FE::new(4) +
            // beta * ry is 8
            test_circuit.y_mid()[0].evaluate(r5) * FE::new(8))
                .representative()
            )
        )
    }

    #[test]
    fn verification_key_gvks_has_length_6_for_test_circuit() {
        let (_, vk): (EvaluationKey<FE>, VerificationKey<FE>) =
            setup(&new_test_qap(), &identity_toxic_waste());
        assert_eq!(vk.gv_ks.len(), 6);
        assert_eq!(vk.gw_ks.len(), 6);
        assert_eq!(vk.gy_ks.len(), 6);
    }
}
