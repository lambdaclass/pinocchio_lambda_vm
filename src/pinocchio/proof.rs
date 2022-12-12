use crate::math;
use math::field_element::FieldElement;

type GroupElement = FieldElement;


pub struct ToxicWaste {
    r_v: FieldElement,
    r_w: FieldElement,
    s: FieldElement,
    alfa_v: FieldElement,
    alfa_w: FieldElement,
    alfa_y: FieldElement,
    beta: FieldElement,
    gamma: FieldElement
}

impl ToxicWaste {
    pub fn new() -> Self {
        Self {
            r_v: FieldElement::random(),
            r_w: FieldElement::random(),
            s: FieldElement::random(),
            alfa_v: FieldElement::random(),
            alfa_w: FieldElement::random(),
            alfa_y: FieldElement::random(),
            beta: FieldElement::random(),
            gamma: FieldElement::random(),
        }
    }
}


struct EvaluationKey {

}

/*
impl EvaluationKey {
    fn new(toxic_waste: ToxicWaste) -> Self {
        // This will be an elliptic curve point
        let g = GroupElement::generator();

        let g_v : Vec<GroupElement>::new();
        for k in intermediate_variables {
            g_v.push(g.pow(toxic_waste.r_v * Vs[k].evaluate(s));
        }
    }
}
*/
