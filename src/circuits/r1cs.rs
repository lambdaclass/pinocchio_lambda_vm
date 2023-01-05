use crate::math::field_element::FieldElement as FE;

#[derive(Debug, PartialEq, Eq)]
pub enum CreationError {
    VectorsSizeMismatch,
}

/// R1CS representation of an Arithmetic Program
#[derive(Clone, Debug, PartialEq, Eq)]
struct Constraint {
    a: Vec<FE>,
    b: Vec<FE>,
    c: Vec<FE>,
}
#[derive(Clone, Debug, PartialEq, Eq)]
/// R1CS represented as a vector of constraints/gates
/// Noticing joining all the first vectors of the constraints results
/// in the first Matrix of the R1CS
/// joining the second vectors results on the second matrix, and so on
struct R1CS {
    constraints: Vec<Constraint>,
}

impl R1CS {
    #[allow(dead_code)]
    pub fn new(constraints: Vec<Constraint>) -> Self {
        Self { constraints }
    }

    #[allow(dead_code)]
    pub fn verify_solution(self, s: &[FE]) -> bool {
        for constraint in self.constraints {
            if !constraint.verify_solution(s) {
                return false;
            }
        }
        true
    }
}

impl Constraint {
    /// Creates a new constraint for a,b,c vectors
    #[allow(dead_code)]
    pub fn new(a: Vec<FE>, b: Vec<FE>, c: Vec<FE>) -> Result<Self, CreationError> {
        if a.len() != b.len() || a.len() != c.len() || b.len() != c.len() {
            Err(CreationError::VectorsSizeMismatch)
        } else {
            Ok(Self { a, b, c })
        }
    }

    #[allow(dead_code)]
    pub fn verify_solution(self, s: &[FE]) -> bool {
        multiply_vectors(&self.a, s) * multiply_vectors(&self.b, s) == multiply_vectors(&self.c, s)
    }
}

pub fn multiply_vectors(v1: &[FE], v2: &[FE]) -> FE {
    let mut acc = FE::new(0);
    for i in 0..v1.len() {
        acc += v1[i] * v2[i]
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mul_vectors_2_2_3_3_equals_12() {
        let v1 = &[FE::new(2), FE::new(2)];
        let v2 = &[FE::new(3), FE::new(3)];

        assert_eq!(multiply_vectors(v1, v2), FE::new(12));
    }

    #[test]
    fn mul_vectors_3_5_equals_15() {
        let v1 = &[FE::new(3)];
        let v2 = &[FE::new(5)];

        assert_eq!(multiply_vectors(v1, v2), FE::new(15));
    }

    #[test]
    fn verify_solution_with_test_circuit_c5_constraints() {
        assert!(test_second_constraint().verify_solution(&test_solution()));
    }

    #[test]
    fn verify_solution_with_test_circuit_c6_constraints() {
        assert!(test_second_constraint().verify_solution(&test_solution()));
    }

    #[test]
    fn verify_bad_solution_with_test_circuit_c5_constraints() {
        let solution = vec![
            FE::new(0),
            FE::new(0),
            FE::new(0),
            FE::new(4),
            FE::new(5),
            FE::new(0),
            FE::new(0),
        ];
        assert!(!test_first_constraint().verify_solution(&solution));
    }

    #[test]
    fn verify_bad_solution_with_test_circuit_c6_constraints() {
        let solution = vec![
            FE::new(0),
            FE::new(2),
            FE::new(1),
            FE::new(4),
            FE::new(5),
            FE::new(2),
            FE::new(2),
        ];
        assert!(!test_second_constraint().verify_solution(&solution));
    }

    #[test]
    fn verify_solution_with_test_r1cs() {
        assert!(test_r1cs().verify_solution(&test_solution()))
    }

    #[test]
    fn verify_bad_solution_with_test_r1cs() {
        let solution = vec![
            FE::new(0),
            FE::new(2),
            FE::new(1),
            FE::new(4),
            FE::new(5),
            FE::new(2),
            FE::new(2),
        ];

        assert!(!test_r1cs().verify_solution(&solution))
    }

    #[test]
    fn verify_bad_solution_because_of_second_constraint_with_test_r1cs() {
        let solution = vec![
            FE::new(0),  // c0
            FE::new(2),  // c1
            FE::new(1),  // c2
            FE::new(5),  // c3
            FE::new(10), // c4
            FE::new(50), // c5 = c4 * c3
            FE::new(2),  // c6 != c5 * (c1+c2), so this should fail
        ];
        assert!(!test_r1cs().verify_solution(&solution))
    }

    #[test]
    fn verify_bad_solution_because_of_first_constraint_with_test_r1cs() {
        let solution = vec![
            FE::new(0),  // c0
            FE::new(1),  // c1
            FE::new(1),  // c2
            FE::new(5),  // c3
            FE::new(10), // c4
            FE::new(10), // c5 != c4 * c3
            FE::new(20), // c6 = c5 * (c1+c2), so this should fail
        ];
        assert!(!test_r1cs().verify_solution(&solution))
    }

    fn test_solution() -> Vec<FE> {
        vec![
            FE::new(0),
            FE::new(1),
            FE::new(2),
            FE::new(3),
            FE::new(4),
            FE::new(12),
            FE::new(36),
        ]
    }

    fn test_r1cs() -> R1CS {
        let constraints = vec![test_first_constraint(), test_second_constraint()];

        R1CS::new(constraints)
    }
    fn test_first_constraint() -> Constraint {
        Constraint {
            a: vec![
                FE::new(0),
                FE::new(0),
                FE::new(0),
                FE::new(1),
                FE::new(0),
                FE::new(0),
                FE::new(0),
            ],
            b: vec![
                FE::new(0),
                FE::new(0),
                FE::new(0),
                FE::new(0),
                FE::new(1),
                FE::new(0),
                FE::new(0),
            ],
            c: vec![
                FE::new(0),
                FE::new(0),
                FE::new(0),
                FE::new(0),
                FE::new(0),
                FE::new(1),
                FE::new(0),
            ],
        }
    }

    fn test_second_constraint() -> Constraint {
        Constraint {
            a: vec![
                FE::new(0),
                FE::new(1),
                FE::new(1),
                FE::new(0),
                FE::new(0),
                FE::new(0),
                FE::new(0),
            ],
            b: vec![
                FE::new(1),
                FE::new(0),
                FE::new(0),
                FE::new(0),
                FE::new(0),
                FE::new(1),
                FE::new(0),
            ],
            c: vec![
                FE::new(0),
                FE::new(0),
                FE::new(0),
                FE::new(0),
                FE::new(0),
                FE::new(0),
                FE::new(1),
            ],
        }
    }
}
