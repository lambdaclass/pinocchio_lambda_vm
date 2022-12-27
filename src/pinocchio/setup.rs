use crate::circuits::qap::Qap;
use crate::math;
use math::field_element::FieldElement as FE;
use math::group::Group;
type GroupType = FE;

#[allow(dead_code)]
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

pub fn setup(qap: Qap, tw: ToxicWaste) -> EvaluationKey {
    let (vs_mid, ws_mid, ys_mid) = (qap.v_mid(), qap.w_mid(), qap.y_mid());

    let s = tw.s;
    let alpha_v = tw.alpha_v;
    let alpha_w = tw.alpha_w;
    let alpha_y = tw.alpha_y;
    let beta = tw.beta;
    let rv = tw.rv;
    let rw = tw.rw;
    let ry = tw.ry();

    let g = FE::generator();

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
