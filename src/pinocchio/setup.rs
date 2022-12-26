use crate::circuits::qap::Qap;
use crate::math;
use math::field_element::FieldElement as FE;
use math::group::Group;
use math::polynomial::Polynomial;
type GroupType = FE;

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

pub fn setup(qap: Qap) -> EvaluationKey {
    let (vs, ws, ys) = (qap.v, qap.w, qap.y);
    
    let s = FE::random();
    let alpha_v = FE::random();
    let alpha_w = FE::random();
    let alpha_y = FE::random();
    let beta = FE::random();
    let rv = FE::random();
    let rw = FE::random();
    let ry = rv * rw;
    let g = FE::generator();
    let ivs = vec![5_usize];
    let degree = 10;

    EvaluationKey {
        gv_ks: ivs
            .iter()
            .map(|&k| g.mul_by_scalar(rv * vs[k].evaluate(s)))
            .collect(),
        gw_ks: ivs
            .iter()
            .map(|&k| g.mul_by_scalar(rw * ws[k].evaluate(s)))
            .collect(),
        gy_ks: ivs
            .iter()
            .map(|&k| g.mul_by_scalar(ry * ys[k].evaluate(s)))
            .collect(),
        gv_alphaks: ivs
            .iter()
            .map(|&k| g.mul_by_scalar(rv * alpha_v * vs[k].evaluate(s)))
            .collect(),
        gw_alphaks: ivs
            .iter()
            .map(|&k| g.mul_by_scalar(rw * alpha_w * ws[k].evaluate(s)))
            .collect(),
        gy_alphaks: ivs
            .iter()
            .map(|&k| g.mul_by_scalar(ry * alpha_y * ys[k].evaluate(s)))
            .collect(),
        g_s_i: (1..=degree).map(|i| g.mul_by_scalar(s.pow(i))).collect(),
        g_beta: ivs
            .iter()
            .map(|&k| {
                g.mul_by_scalar(
                    rv * beta * vs[k].evaluate(s)
                        + rw * beta * ws[k].evaluate(s)
                        + ry * beta * ys[k].evaluate(s),
                )
            })
            .collect(),
    }
}
