#![allow(dead_code)]

use na::Matrix3;

pub trait ApproxEqual {
    type Epsilon: Copy;
    fn approx_eq(&self, other: Self, eps: Self::Epsilon) -> bool;
    fn default_epsilon() -> Self::Epsilon;
}

impl ApproxEqual for f64 {
    type Epsilon = f64;

    fn approx_eq(&self, other: Self, eps: Self::Epsilon) -> bool {
        (self - other).abs() < eps
    }

    fn default_epsilon() -> Self::Epsilon {
        1e-8
    }
}

pub fn all_equal_with_epsilon<T: ApproxEqual + Copy>(slice: &[T], eps: T::Epsilon) -> bool {
    match slice.first() {
        None => true,
        Some(&first) => slice.iter().all(|&x| x.approx_eq(first, eps)),
    }
}

pub fn all_equal<T: ApproxEqual + Copy>(slice: &[T]) -> bool {
    match slice.first() {
        None => true,
        Some(&first) => slice
            .iter()
            .all(|&x| x.approx_eq(first, T::default_epsilon())),
    }
}

pub fn symmetrize(a_matrix: &Matrix3<f64>) -> Matrix3<f64> {
    (a_matrix + a_matrix.transpose()) * 0.5
}

pub fn mat_exp_taylor(x: &Matrix3<f64>) -> Matrix3<f64> {
    let x2 = x * x;
    let x3 = x2 * x;
    let x4 = x3 * x;
    Matrix3::identity() + x + x2 * (1.0 / 2.0) + x3 * (1.0 / 6.0) + x4 * (1.0 / 24.0)
}
