use crate::math::cyclic_group::CyclicBilinearGroup;
use crate::math::elliptic_curve::EllipticCurveElement;
use crate::math::field_element::FieldElement;

type FE = FieldElement<5>;
pub type CyclicBilinearGroupType = EllipticCurveElement;

/// Calculates msm for C and hidings
/// if either array is empty, returns zero
pub fn msm<T>(c: &[FE], hidings: &[T]) -> T
where
    T: CyclicBilinearGroup,
{
    c.iter()
        .zip(hidings.iter())
        .map(|(&c, h)| h.operate_with_self(c.representative()))
        .reduce(|acc, x| acc.operate_with(&x))
        .unwrap_or_else(T::neutral_element)
}

#[cfg(test)]
mod tests {
    use super::*;

    //MSM tests require the CyclicBilinearGroupType to be a FieldElement
    #[test]
    fn msm_11_is_1() {
        let c = [FE::new(1)];
        let hiding = [CyclicBilinearGroupType::generator()];
        assert_eq!(msm(&c, &hiding), CyclicBilinearGroupType::generator());
    }

    #[test]
    fn msm_23_is_6() {
        let c = [FE::new(3)];
        let hiding = [FE::new(2)];
        assert_eq!(msm(&c, &hiding), FE::new(6));
    }

    #[test]
    fn msm_with_c_2_3_hiding_3_4_is_18() {
        let c = [FE::new(2), FE::new(3)];
        let hiding = [FE::new(3), FE::new(4)];
        assert_eq!(msm(&c, &hiding), FE::new(18));
    }

    #[test]
    fn msm_with_empty_c_is_none() {
        let c = [];
        let hiding = [FE::new(3), FE::new(4)];
        assert_eq!(msm(&c, &hiding), FE::new(0));
    }

    #[test]
    fn msm_with_emtpy_hiding_is_none() {
        let c = [FE::new(0)];
        let hiding: Vec<FE> = Vec::new();
        assert_eq!(msm(&c, &hiding), FE::new(0));
    }

    #[test]
    fn msm_with_empty_arguments_is_none() {
        let c = [];
        let hiding: Vec<FE> = Vec::new();
        assert_eq!(msm(&c, &hiding), FE::new(0));
    }
}
