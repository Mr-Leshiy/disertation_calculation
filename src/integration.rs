use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use std::f64::consts::PI;

// calculate a definite integral of the function on the provided interval
pub fn definite_integral<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    mut n: u32,
    eps: f64,
    f: &F,
) -> f64 {
    assert!(a < b);

    let mut h = (b - a) / n as f64;
    let mut prev_result;

    let mut result = (0..n)
        .into_par_iter()
        .map(|i| {
            let x = a + (i as f64) * h;
            f(x)
        })
        .sum();

    result *= h;

    loop {
        n *= 2;
        h /= 2_f64;
        prev_result = result;

        result = (0..n)
            .into_par_iter()
            .map(|i| {
                let x = a + (i as f64) * h;
                f(x)
            })
            .sum();
        result *= h;

        if f64::abs(result - prev_result) < eps {
            break;
        }
    }

    result
}

#[allow(dead_code)]
pub fn definite_integral_limit<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    n: u32,
    f: &F,
) -> f64 {
    assert!(a < b);

    let h = (b - a) / n as f64;
    let mut result = (0..n)
        .into_par_iter()
        .map(|i| {
            let x = a + (i as f64) * h;
            f(x)
        })
        .sum();
    result *= h;

    result
}

pub fn sqrt_gauss_integral<F: Fn(f64) -> f64 + Send + Sync>(mut n: u32, eps: f64, f: &F) -> f64 {
    let mut prev_result;

    let mut result = (1..n)
        .into_par_iter()
        .map(|i| {
            let x = f64::cos(PI * (2_f64 * i as f64 - 1_f64) / (2_f64 * n as f64));
            f(x)
        })
        .sum::<f64>()
        * PI
        / n as f64;

    loop {
        n *= 2;
        prev_result = result;

        result = (1..n)
            .into_par_iter()
            .map(|i| {
                let x = f64::cos(PI * (2_f64 * i as f64 - 1_f64) / (2_f64 * n as f64));
                f(x)
            })
            .sum::<f64>()
            * PI
            / n as f64;

        if f64::abs(result - prev_result) < eps {
            break;
        }
    }
    result
}

pub fn sqrt_gauss_integral_finit<F: Fn(f64) -> f64 + Send + Sync>(n: u32, f: &F) -> f64 {
    let result = (1..n)
        .into_par_iter()
        .map(|i| {
            let x = f64::cos(PI * (2_f64 * i as f64 - 1_f64) / (2_f64 * n as f64));
            f(x)
        })
        .sum::<f64>()
        * PI
        / n as f64;
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polynomials::chebyshev;

    #[test]
    fn test_integration() {
        let eps = 0.00001;

        assert!(
            f64::abs(definite_integral(0_f64, 2_f64, 10, eps, &|x| x) - 2_f64) < eps,
            "result, {}",
            definite_integral(0_f64, 2_f64, 10, eps, &|x| x)
        );
        assert!(
            f64::abs(definite_integral(0_f64, 2_f64, 10, eps, &|x| x * x) - 8_f64 / 3_f64) < eps,
            "result, {}",
            definite_integral(0_f64, 2_f64, 10, eps, &|x| x * x)
        );
        assert!(
            f64::abs(
                definite_integral(0_f64, 2_f64, 10, eps, &|x| f64::sin(x) * 5_f64)
                    - (-5_f64 * f64::cos(2_f64) + 5_f64)
            ) < eps,
            "result, {}",
            definite_integral(0_f64, 2_f64, 10, eps, &|x| f64::sin(x) * 5_f64)
        );
    }

    #[test]
    fn spectral_relation_test() {
        let x = 2.0;
        let x = f64::cos(PI / 2.0 * x);
        let m = 2;
        let eps = 0.001;

        let f = |y| {
            let cheb = chebyshev(y, 2 * m + 1);
            let log = f64::ln(f64::abs((x + y) / (x - y)));
            cheb * log
        };

        let right = sqrt_gauss_integral(10, eps, &f) / std::f64::consts::PI / 2_f64;
        let left = chebyshev(x, 2 * m + 1) / (2_f64 * m as f64 + 1_f64);
        println!("right: {right}, left: {left}");
    }
}
