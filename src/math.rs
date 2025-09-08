#![allow(dead_code)]

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
        Some(&first) => slice.iter().all(|&x| x.approx_eq(first, T::default_epsilon()))
        
    }
}