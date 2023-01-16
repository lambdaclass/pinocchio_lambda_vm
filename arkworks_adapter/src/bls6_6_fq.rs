use ark_ff::{biginteger::BigInteger64 as BigInteger, fields::*};

pub type Fq = Fp64<FqParameters>;

pub struct FqParameters;

impl Fp64Parameters for FqParameters {}

impl FftParameters for FqParameters {
    type BigInt = BigInteger;

    // N = 2^s * t, t = 3 (odd integer), N = 4
    // -> s = 3
    const TWO_ADICITY: u32 = 3;

    #[rustfmt::skip]
    const TWO_ADIC_ROOT_OF_UNITY: BigInteger = BigInteger([
        3
    ]);
}

impl FpParameters for FqParameters {
    const MODULUS: BigInteger = BigInteger([5]);

    const MODULUS_BITS: u32 = 3;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const REPR_SHAVE_BITS: u32 = 0;

    const R: BigInteger = BigInteger([1]);

    const R2: BigInteger = BigInteger([1]);

    const INV: u64 = 3689348814741910323;

    /// GENERATOR = 5
    /// Encoded in Montgomery form, so the value here is
    #[rustfmt::skip]
    const GENERATOR: BigInteger = BigInteger([3]);

    #[rustfmt::skip]
    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([1]);

    // This is used for Tonelly Shanks, an algorithm to solve a square root.

    // 2^s * t=5-1 -> 2^s * 3 = 4
    // -> s = 
    #[rustfmt::skip]
    const T: BigInteger = BigInteger([3]);

    #[rustfmt::skip]
    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([2]);
}

#[allow(dead_code)]
pub const FQ_ONE: Fq = Fq::new(FqParameters::R);
#[allow(dead_code)]
pub const FQ_ZERO: Fq = Fq::new(BigInteger([0]));
