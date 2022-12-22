struct qap {
    v: Vec<Polynomial>,
    w: Vec<Polynomial>,
    y: Vec<Polynomial>,
    target: Polynomial
}

#[derive(Debug, PartialEq, Eq)]
pub enum CreationError {
    VectorsSizeMismatch
}

impl qap {
    pub fn new(
        v: Vec<Polynomial>,
        w: Vec<Polynomial>,
        y: Vec<Polynomial>,
        target: Polynomial) 
    -> Result<Self, PolynomialSizeMismatch> {
        if v.len() != w.len() | v.len() != y.len() | w.len() != y.len() {
            PolynomialSizeMismatch
        } else {
            qap {
                v,
                w,
                y,
                target
            }
        }
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn qap_with_different_amount_of_polynomials_should_error() -> Polynomial {

        let v = vec![Polynomial::new(vec![FE::new(1).unwrap(),FE::new(2).unwrap()])];
        let u = v.clone();
        let w = vec![Polynomial::new(vec![FE::new(2).unwrap()]); 
        let t = Polynomial::new(vec![FE::new(3).unwrap()]);
        assert_eq!(qap.new(v,u,w,t), Err(_));
    }
}
