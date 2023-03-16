use rayon::prelude::{IntoParallelIterator, ParallelIterator};

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
        .into_iter()
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

#[cfg(test)]
mod tests {
    use super::*;

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
}
