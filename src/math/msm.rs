use crate::math::field_element::FieldElement as FE;

use crate::math::cyclic_group::CyclicGroup;

pub type CyclicGroupType = FE;

/// Calculates msm for C and hidings
/// if either array is empty, returns zero
pub fn msm(c: &[FE], hidings: &[CyclicGroupType]) -> CyclicGroupType {
    c.iter()
        .zip(hidings.iter())
        .map(|(&c, &h)| h.operate_with_self(c))
        .reduce(|acc, x| acc + x)
        .unwrap_or_else(CyclicGroupType::zero)
}

#[cfg(test)]
mod tests {
    use super::*;

    //MSM tests require the CyclicGroupType to be a FieldElement
    #[test]
    fn msm_11_is_1() {
        let c = [FE::one()];
        let hiding = [FE::one()];
        assert_eq!(msm(&c, &hiding), FE::one());
    }

    #[test]
    fn msm_23_is_6() {
        let c = [FE::new(3).unwrap()];
        let hiding = [FE::new(2).unwrap()];
        assert_eq!(msm(&c, &hiding), FE::new(6).unwrap());
    }

    #[test]
    fn msm_with_c_2_3_hiding_3_4_is_18() {
        let c = [FE::new(2).unwrap(), FE::new(3).unwrap()];
        let hiding = [FE::new(3).unwrap(), FE::new(4).unwrap()];
        assert_eq!(msm(&c, &hiding), FE::new(18).unwrap());
    }

    #[test]
    fn msm_with_empty_c_is_none() {
        let c = [];
        let hiding = [FE::new(3).unwrap(), FE::new(4).unwrap()];
        assert_eq!(msm(&c, &hiding), FE::zero());
    }

    #[test]
    fn msm_with_emtpy_hiding_is_none() {
        let c = [FE::zero()];
        let hiding = [];
        assert_eq!(msm(&c, &hiding), FE::zero());
    }

    #[test]
    fn msm_with_empty_arguments_is_none() {
        let c = [];
        let hiding = [];
        assert_eq!(msm(&c, &hiding), FE::zero());
    }
}
