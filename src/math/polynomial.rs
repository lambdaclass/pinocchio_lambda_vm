use super::field_element::FieldElement as FE;

/// Evaluate "e" on the Lagrange Polynomial interpolated on the field elements (xs, ys)
pub fn interpolate(xs: &[FE], ys: &[FE], e: FE) -> FE {
    let mut result: FE = FE::new(0).unwrap();

    for (i, y) in ys.iter().enumerate() {
        let mut y_term: FE = *y;
        for (j, x) in xs.iter().enumerate() {
            if i != j {
                y_term = y_term * ((e - *x) / (xs[i] - *x));
            }
        }
        result = result + y_term;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interpolate_x_0_2_y_3_4() {
        let x = interpolate(
            &[FE::new(0).unwrap(), FE::new(2).unwrap()],
            &[FE::new(3).unwrap(), FE::new(4).unwrap()],
            FE::new(2).unwrap(),
        );
        assert_eq!(FE::new(4).unwrap(), x);
    }

    #[test]
    fn interpolate_x_2_5_7_y_10_19_43() {
        let xs = &[
            FE::new(2).unwrap(),
            FE::new(5).unwrap(),
            FE::new(7).unwrap(),
        ];
        let ys = &[
            FE::new(10).unwrap(),
            FE::new(19).unwrap(),
            FE::new(43).unwrap(),
        ];

        assert_eq!(
            FE::new(10).unwrap(),
            interpolate(xs, ys, FE::new(2).unwrap())
        );
        assert_eq!(
            FE::new(19).unwrap(),
            interpolate(xs, ys, FE::new(5).unwrap())
        );
        assert_eq!(
            FE::new(43).unwrap(),
            interpolate(xs, ys, FE::new(7).unwrap())
        );
    }
}
