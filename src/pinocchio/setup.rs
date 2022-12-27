use crate::circuits::qap::Qap;
use crate::math;
use math::field_element::FieldElement as FE;
use math::group::Group;
type GroupType = FE;

#[allow(dead_code)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EvaluationKey {
    gv_ks: Vec<GroupType>,
    gw_ks: Vec<GroupType>,
    gy_ks: Vec<GroupType>,
    gv_alphaks: Vec<GroupType>,
    gw_alphaks: Vec<GroupType>,
    gy_alphaks: Vec<GroupType>,
    g_s_i: Vec<GroupType>,
    g_beta: Vec<GroupType>,
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
}

impl ToxicWaste {
    pub fn ry(self) -> FE {
        self.rv * self.rw
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
        }
    }
}

pub fn setup(qap: Qap, toxic_waste: ToxicWaste) -> EvaluationKey {
    let (vs_mid, ws_mid, ys_mid) = (qap.v_mid(), qap.w_mid(), qap.y_mid());

    let s = toxic_waste.s;
    let alpha_v = toxic_waste.alpha_v;
    let alpha_w = toxic_waste.alpha_w;
    let alpha_y = toxic_waste.alpha_y;
    let beta = toxic_waste.beta;
    let rv = toxic_waste.rv;
    let rw = toxic_waste.rw;
    let ry = toxic_waste.ry();

    let g = GroupType::generator();

    let degree = qap.target.degree();

    let mut gv_ks: Vec<GroupType> = Vec::with_capacity(vs_mid.len());
    let mut gw_ks: Vec<GroupType> = Vec::with_capacity(vs_mid.len());
    let mut gy_ks: Vec<GroupType> = Vec::with_capacity(vs_mid.len());
    let mut gv_alphaks: Vec<GroupType> = Vec::with_capacity(vs_mid.len());
    let mut gw_alphaks: Vec<GroupType> = Vec::with_capacity(vs_mid.len());
    let mut gy_alphaks: Vec<GroupType> = Vec::with_capacity(vs_mid.len());
    let mut g_beta: Vec<GroupType> = Vec::with_capacity(vs_mid.len());
    // g_s_i is the only paramater to depend on the degree of the qap
    // This is an upper bound, it could be smaller
    let mut g_s_i: Vec<GroupType> = Vec::with_capacity(degree);

    // Set evaluation keys for each of their respective k mid element
    for k in 0..vs_mid.len() {
        gv_ks.push(g.mul_by_scalar(rv * vs_mid[k].evaluate(s)));
        gw_ks.push(g.mul_by_scalar(rw * ws_mid[k].evaluate(s)));
        gy_ks.push(g.mul_by_scalar(ry * ys_mid[k].evaluate(s)));
        gv_alphaks.push(g.mul_by_scalar(ry * alpha_v * vs_mid[k].evaluate(s)));
        gw_alphaks.push(g.mul_by_scalar(rw * alpha_w * ws_mid[k].evaluate(s)));
        gy_alphaks.push(g.mul_by_scalar(ry * alpha_y * ys_mid[k].evaluate(s)));
        g_beta.push(
            rv * beta * vs_mid[k].evaluate(s)
                + rw * beta * ws_mid[k].evaluate(s)
                + ry * beta * ys_mid[k].evaluate(s),
        )
    }

    for i in 0..qap.target.degree() {
        // This unwrap would only fail in an OS
        // with 256 bits pointer, which doesn't exist
        g_s_i.push(g.mul_by_scalar(s.pow(i.try_into().unwrap())));
    }

    EvaluationKey {
        gv_ks,
        gw_ks,
        gy_ks,
        gv_alphaks,
        gw_alphaks,
        gy_alphaks,
        g_s_i,
        g_beta,
    }
}

#[cfg(test)]
mod tests {
    use super::{setup, ToxicWaste};
    use crate::{circuits::qap::Qap, math};
    use math::field_element::FieldElement as FE;

    fn identity_toxic_waste() -> ToxicWaste {
        ToxicWaste {
            s: FE::one(),
            alpha_v: FE::one(),
            alpha_w: FE::one(),
            alpha_y: FE::one(),
            beta: FE::one(),
            rv: FE::one(),
            rw: FE::one(),
        }
    }

    #[test]
    fn evaluation_keys_size_for_test_circuit_is_1_for_each_key() {
        let eval_key = setup(Qap::new_test_circuit(), identity_toxic_waste());
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
        let r5 = FE::new(2).unwrap();
        let tw = ToxicWaste {
            s: FE::new(2).unwrap(),
            alpha_v: FE::new(2).unwrap(),
            alpha_w: FE::new(2).unwrap(),
            alpha_y: FE::new(2).unwrap(),
            beta: FE::new(2).unwrap(),
            rv: FE::new(2).unwrap(),
            rw: FE::new(2).unwrap(),
        };

        let test_circuit = Qap::new_test_circuit();

        let eval_key = setup(test_circuit.clone(), tw);

        // These keys should be the same evaluation * rv, which is two
        assert_eq!(
            eval_key.gv_ks[0],
            test_circuit.v_mid()[0].evaluate(r5) * FE::new(2).unwrap()
        );
        assert_eq!(
            eval_key.gw_ks[0],
            test_circuit.w_mid()[0].evaluate(r5) * FE::new(2).unwrap()
        );
        // These keys should be the same evaluation * ys, which is two
        // Since the whole thing is 0
        assert_eq!(
            eval_key.gy_ks[0],
            test_circuit.y_mid()[0].evaluate(r5) * FE::new(4).unwrap()
        );

        // alpha * rv and alpha * rw is 4
        assert_eq!(
            eval_key.gv_alphaks[0],
            test_circuit.v_mid()[0].evaluate(r5) * FE::new(4).unwrap()
        );
        assert_eq!(
            eval_key.gv_alphaks[0],
            test_circuit.v_mid()[0].evaluate(r5) * FE::new(4).unwrap()
        );
        // alpha * ry and alpha * rw is 8
        assert_eq!(
            eval_key.gv_alphaks[0],
            test_circuit.v_mid()[0].evaluate(r5) * FE::new(8).unwrap()
        );

        assert_eq!(
            eval_key.g_beta[0],
            // beta * rv is 4
            test_circuit.v_mid()[0].evaluate(r5) * FE::new(4).unwrap() +
            test_circuit.w_mid()[0].evaluate(r5) * FE::new(4).unwrap() +
            // beta * ry is 8
            test_circuit.y_mid()[0].evaluate(r5) * FE::new(8).unwrap()
        )
    }
}
