pub fn chebyshev(x: f64, n: usize) -> f64 {
    if n == 0 {
        return 1_f64;
    }
    if n == 1 {
        return x;
    }

    let mut t_n_1 = 1_f64;
    let mut t_n_2 = x;
    let mut t_n_3 = 2_f64 * x * t_n_2 - t_n_1;
    for _ in 2..n {
        t_n_1 = t_n_2;
        t_n_2 = t_n_3;
        t_n_3 = 2_f64 * x * t_n_2 - t_n_1;
    }

    t_n_3
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chebyshev_test() {
        for x in 0..100 {
            let x = x as f64;
            assert_eq!(chebyshev(x, 0), 1_f64);
            assert_eq!(chebyshev(x, 1), x);
            assert_eq!(chebyshev(x, 2), 2_f64 * x * x - 1_f64);
            assert_eq!(chebyshev(x, 3), 4_f64 * x * x * x - 3_f64 * x);
            assert_eq!(
                chebyshev(x, 4),
                8_f64 * x * x * x * x - 8_f64 * x * x + 1_f64
            );
            assert_eq!(
                chebyshev(x, 5),
                16_f64 * x * x * x * x * x - 20_f64 * x * x * x + 5_f64 * x
            );
            assert_eq!(
                chebyshev(x, 6),
                32_f64 * x * x * x * x * x * x - 48_f64 * x * x * x * x + 18_f64 * x * x - 1_f64
            );
            assert_eq!(
                chebyshev(x, 7),
                64_f64 * x * x * x * x * x * x * x - 112_f64 * x * x * x * x * x
                    + 56_f64 * x * x * x
                    - 7_f64 * x
            );
            assert_eq!(
                chebyshev(x, 8),
                128_f64 * x * x * x * x * x * x * x * x - 256_f64 * x * x * x * x * x * x
                    + 160_f64 * x * x * x * x
                    - 32_f64 * x * x
                    + 1_f64
            );
        }
    }
}
