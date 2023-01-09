/// Cyclic group implementation for Pinocchio Implementation
/// inverse function is not implemented since it's not needed for it
pub trait CyclicGroup {
    fn generator() -> Self;
    fn neutral_element() -> Self;
    /// Repeats the group operation with self `times` times
    /// Operation can be add or mul depending on the GroupType
    fn operate_with_self(&self, times: u128) -> Self;
    /// Operation can be add or mul depending on the GroupType
    fn operate_with(&self, other: &Self) -> Self;
    fn pairing(&self, other: &Self) -> Self;
}
