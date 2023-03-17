use crate::{
    integration::definite_integral,
    polynomials::chebyshev,
    utils::{g, lambda, mu_0},
};
use std::f64::consts::PI;

fn h_m(a: f64, m: usize, eps: f64) -> f64 {
    let f = |x| {
        let a1 = chebyshev(f64::cos(x * PI / (2_f64 * a)), 2 * m + 1);
        let a2 = f64::sqrt(1_f64 - x * x);
        a1 / a2
    };
    definite_integral(-1_f64, 1_f64, 100, eps, &f)
}

fn f_m<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    m: usize,
    eps: f64,
) -> f64 {
    let a4 = a * (f64::cosh(b * PI / a) - 1_f64);
    // let pn = definite_integral(0_f64, a, 100, eps, &|x| {
    //     load_function(x) * f64::cos(alpha * x)
    // });
    let f = |x| {
        let a1 = chebyshev(f64::cos(x * PI / (2_f64 * a)), 2 * m + 1);
        let a2 = f64::sqrt(1_f64 - x * x);
        a1 / a2
    };
    definite_integral(-1_f64, 1_f64, 100, eps, &f)
}

fn coefficients(
    b: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
) -> (
    (f64, f64, f64, f64, f64, f64, f64, f64),
    (f64, f64, f64, f64, f64, f64, f64, f64),
) {
    let coef = (1_f64 + mu_0) * 4_f64 * alpha;
    let e1 = f64::exp(alpha * b);
    let e2 = f64::exp(-alpha * b);

    let x1 = 2_f64 * b * alpha * alpha * mu_0 + 2_f64 * alpha + 2_f64 * alpha * mu_0;
    let x2 = 2_f64 * b * alpha * alpha * mu_0 - 2_f64 * alpha;
    let x3 = -2_f64 * b * alpha * alpha * mu_0 + 2_f64 * alpha + 2_f64 * alpha * mu_0;
    let x4 = 2_f64 * b * alpha * alpha * mu_0 + 2_f64 * alpha;
    let x5 =
        -2_f64 * g * b * alpha * alpha * mu_0 - 2_f64 * g * alpha * mu_0 + 2_f64 * alpha * lambda;
    let x6 = -2_f64 * g * b * alpha * alpha * mu_0 + 4_f64 * g * alpha + 2_f64 * alpha * lambda;
    let x7 =
        -2_f64 * g * b * alpha * alpha * mu_0 + 2_f64 * g * alpha * mu_0 - 2_f64 * alpha * lambda;
    let x8 = 2_f64 * g * b * alpha * alpha * mu_0 + 4_f64 * g * alpha + 2_f64 * alpha * lambda;
    let x9 = 2_f64 * alpha + 2_f64 * alpha * mu_0;
    let x10 = -2_f64 * alpha;
    let x11 = 2_f64 + mu_0;
    let x12 = -2_f64 - mu_0;

    let y1 =
        e1 * (x1 * x10 * x12) / (x9 * x11) + e1 * (x1 * x10) / x9 - e1 * (x2 * x12) / x11 + e2 * x4;
    let y2 = e2 * x3 - e1 * x1;
    let y3 =
        e1 * (x5 * x10 * x12) / (x9 * x11) + e1 * (x5 * x10) / x9 - e1 * (x6 * x12) / x11 + e2 * x8;
    let y4 = e2 * x7 - e1 * x5;

    let d_1_0 = -coef * (y4 * x10 * x12 + y4 * x10 * x11 + y3 * x9 * x11)
        / (x9 * x11 * (y2 * y3 - y1 * y4));
    let d_2_0 = -coef * (y2 * x10 * x12 + y2 * x10 * x11 + y1 * x9 * x11)
        / (x9 * x11 * (y4 * y1 - y3 * y2));
    let d_3_0 = coef * (x12 * y4) / (x11 * (y2 * y3 - y1 * y4));
    let d_4_0 = coef * (x12 * y2) / (x11 * (y4 * y1 - y3 * y2));
    let f_1_0 = coef * y3 / (y2 * y3 - y1 * y4);
    let f_2_0 = coef * y1 / (y4 * y1 - y3 * y2);
    let f_3_0 = -coef * y4 / (y2 * y3 - y1 * y4);
    let f_4_0 = -coef * y2 / (y4 * y1 - y3 * y2);

    let d_1_1 = coef
        * e1
        * ((y2 * (x1 * y4 - x5 * y2) * (x10 * x12 + x10 * x11)
            + x11 * x9 * (x1 * y4 - x5 * y2) * y1
            + x11 * x9 * x1 * (y3 * y2 - y1 * y4))
            / (x11 * x9 * x9 * y2 * (y3 * y2 - y1 * y4)))
        + coef / x9;
    let d_3_1 = -coef * e1 * (x12 * (x1 * y4 - x5 * y2) / (x11 * x9 * (y3 * y2 - y1 * y4)));
    let f_1_1 = -coef
        * e1
        * ((x1 * (y3 * y2 - y1 * y4) + (x1 * y4 - x5 * y2) * y1) / (x9 * y2 * (y3 * y2 - y1 * y4)));
    let f_3_1 = coef * e1 * (x1 * y4 - x5 * y2) / (x9 * (y3 * y2 - y1 * y4));

    let d_2_1 = coef
        * e1
        * (((x5 * x10 - x1 * x9) * y2 - (x1 * x10 - x2 * x9) * y4) * (x10 * x11 + x10 * x12)
            / (x9 * x9 * x11 * x11 * (y2 * y3 - y1 * y4))
            - ((x1 * x10 - x2 * x9) * y3 - (x5 * x10 - x6 * x9) * y1)
                / (x9 * x11 * (y2 * y3 - y1 * y4)))
        - coef * x10 / (x9 * x11);
    let d_4_1 = -coef * e1 * ((x5 * x10 - x6 * x9) * y2 - (x1 * x10 - x2 * x9) * y4) * x12
        / (x9 * x11 * x11 * (y2 * y3 - y1 * y4))
        + coef / x11;
    let f_2_1 = coef * e1 * ((x1 * x10 - x2 * x9) * y3 - (x5 * x10 - x6 * x9) * y1)
        / (x9 * x11 * (y2 * y3 - y1 * y4));
    let f_4_1 = coef * e1 * ((x5 * x10 - x6 * x9) * y2 - (x1 * x10 - x2 * x9) * y4)
        / (x9 * x11 * (y2 * y3 - y1 * y4));

    (
        (d_1_0, d_2_0, d_3_0, d_4_0, f_1_0, f_2_0, f_3_0, f_4_0),
        (d_1_1, d_2_1, d_3_1, d_4_1, f_1_1, f_2_1, f_3_1, f_4_1),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn check() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;
        let eps = 0.0001;

        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);

        let n = 1000000000_f64;
        let alpha = PI * n / a;
        let e2 = f64::exp(-alpha * b);

        let x1 = 2_f64 * b * mu_0 + 2_f64 / alpha + 2_f64 * mu_0 / alpha;
        let x2 = 2_f64 * b * mu_0 - 2_f64 / alpha;
        let x3 = -2_f64 * b * mu_0 + 2_f64 / alpha + 2_f64 * mu_0 / alpha;
        let x4 = 2_f64 * b * mu_0 + 2_f64 / alpha;
        let x5 = -2_f64 * g * b * mu_0 - 2_f64 * g * mu_0 / alpha + 2_f64 * lambda / alpha;
        let x6 = -2_f64 * g * b * mu_0 + 4_f64 * g / alpha + 2_f64 * lambda / alpha;
        let x7 = -2_f64 * g * b * mu_0 + 2_f64 * g * mu_0 / alpha - 2_f64 * lambda / alpha;
        let x8 = 2_f64 * g * b * mu_0 + 4_f64 * g / alpha + 2_f64 * lambda / alpha;
        let x9 = 2_f64 + 2_f64 * mu_0;
        let x10 = -2_f64;
        let x11 = 2_f64 + mu_0;
        let x12 = -2_f64 - mu_0;

        let y1 = (x1 * x10 * x12) / (x9 * x11) + (x1 * x10) / x9 - (x2 * x12) / x11 + e2 * x4;
        let y2 = e2 * e2 * x3 - x1;
        let y3 = (x5 * x10 * x12) / (x9 * x11) + (x5 * x10) / x9 - (x6 * x12) / x11 + e2 * e2 * x8;
        let y4 = e2 * e2 * x7 - x5;

        assert!(f64::abs(x1 - 2_f64 * b * mu_0) < eps);
        assert!(f64::abs(x2 - 2_f64 * b * mu_0) < eps);
        assert!(f64::abs(x3 + 2_f64 * b * mu_0) < eps);
        assert!(f64::abs(x4 - 2_f64 * b * mu_0) < eps);
        assert!(f64::abs(x5 + 2_f64 * g * b * mu_0) < eps);
        assert!(f64::abs(x6 + 2_f64 * g * b * mu_0) < eps);
        assert!(f64::abs(x7 + 2_f64 * g * b * mu_0) < eps);
        assert!(f64::abs(x8 - 2_f64 * g * b * mu_0) < eps);
        assert!(f64::abs(x9 - (2_f64 + 2_f64 * mu_0)) < eps);
        assert!(f64::abs(x10 + 2_f64) < eps);
        assert!(f64::abs(x11 - (2_f64 + mu_0)) < eps);
        assert!(f64::abs(x12 - (-2_f64 - mu_0)) < eps);

        assert!(f64::abs(y1 - 2_f64 * b * mu_0) < eps);
        assert!(f64::abs(y2 + 2_f64 * b * mu_0) < eps);
        assert!(f64::abs(y3 + 2_f64 * g * b * mu_0) < eps);
        assert!(f64::abs(y4 - 2_f64 * g * b * mu_0) < eps);

        let coef = (1_f64 + mu_0) * 4_f64 * alpha;

        let d_1_0 = -coef * (y4 * x10 * x12 + y4 * x10 * x11 + y3 * x9 * x11)
            / (x9 * x11 * (y2 * y3 - y1 * y4));
        let d_2_0 = -coef * (y2 * x10 * x12 + y2 * x10 * x11 + y1 * x9 * x11)
            / (x9 * x11 * (y4 * y1 - y3 * y2));
        let d_3_0 = coef * (x12 * y4) / (x11 * (y2 * y3 - y1 * y4));
        let d_4_0 = coef * (x12 * y2) / (x11 * (y4 * y1 - y3 * y2));
        let f_1_0 = coef * y3 / (y2 * y3 - y1 * y4);
        let f_2_0 = coef * y1 / (y4 * y1 - y3 * y2);
        let f_3_0 = -coef * y4 / (y2 * y3 - y1 * y4);
        let f_4_0 = -coef * y2 / (y4 * y1 - y3 * y2);
        println!("d_1_0: {0}, d_2_0: {1}, d_3_0: {2}, d_4_0: {3}, f_1_0: {4}, f_2_0: {5}, f_3_0: {6}, f_4_0: {7}", d_1_0, d_2_0, d_3_0, d_4_0, f_1_0,f_2_0, f_3_0, f_4_0);
    }

    #[test]
    fn coefficients_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;
        let eps = 0.0001;

        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);

        for i in 1..3 {
            let alpha = PI * i as f64 / a;

            let e1 = f64::exp(alpha * b);
            let e2 = f64::exp(-alpha * b);

            let x1 = 2_f64 * b * alpha * alpha * mu_0 + 2_f64 * alpha + 2_f64 * alpha * mu_0;
            let x2 = 2_f64 * b * alpha * alpha * mu_0 - 2_f64 * alpha;
            let x3 = -2_f64 * b * alpha * alpha * mu_0 + 2_f64 * alpha + 2_f64 * alpha * mu_0;
            let x4 = 2_f64 * b * alpha * alpha * mu_0 + 2_f64 * alpha;
            let x5 = -2_f64 * g * b * alpha * alpha * mu_0 - 2_f64 * g * alpha * mu_0
                + 2_f64 * alpha * lambda;
            let x6 =
                -2_f64 * g * b * alpha * alpha * mu_0 + 4_f64 * g * alpha + 2_f64 * alpha * lambda;
            let x7 = -2_f64 * g * b * alpha * alpha * mu_0 + 2_f64 * g * alpha * mu_0
                - 2_f64 * alpha * lambda;
            let x8 =
                2_f64 * g * b * alpha * alpha * mu_0 + 4_f64 * g * alpha + 2_f64 * alpha * lambda;
            let x9 = 2_f64 * alpha + 2_f64 * alpha * mu_0;
            let x10 = -2_f64 * alpha;
            let x11 = 2_f64 + mu_0;
            let x12 = -2_f64 - mu_0;

            let (
                (d_1_0, d_2_0, d_3_0, d_4_0, f_1_0, f_2_0, f_3_0, f_4_0),
                (d_1_1, d_2_1, d_3_1, d_4_1, f_1_1, f_2_1, f_3_1, f_4_1),
            ) = coefficients(b, alpha, mu_0, g, lambda);

            let eq1 = |d1, d2, f1, f2| e1 * x1 * d1 + e1 * x2 * d2 + e2 * x3 * f1 + e2 * x4 * f2;
            let eq2 = |d1, d2, f1, f2| e1 * x5 * d1 + e1 * x6 * d2 + e2 * x7 * f1 + e2 * x8 * f2;
            let eq3 = |d1, d2, f1, f2| x9 * d1 + x10 * d2 + x9 * f1 - x10 * f2;
            let eq4 = |d1, f1| x11 * d1 + x12 * f1;
            //
            assert!(
                f64::abs(eq1(d_1_0, d_3_0, f_1_0, f_3_0) - (1_f64 + mu_0) * 4_f64 * alpha) < eps,
                "result: {}, exp: {}",
                eq1(d_1_0, d_3_0, f_1_0, f_3_0),
                (1_f64 + mu_0) * 4_f64 * alpha
            );
            assert!(
                f64::abs(eq2(d_1_0, d_3_0, f_1_0, f_3_0) - 0_f64) < eps,
                "result: {}, exp: {}",
                eq2(d_1_0, d_3_0, f_1_0, f_3_0),
                0_f64
            );
            assert!(
                f64::abs(eq3(d_1_0, d_3_0, f_1_0, f_3_0) - 0_f64) < eps,
                "result: {}, exp: {}",
                eq3(d_1_0, d_3_0, f_1_0, f_3_0),
                0_f64
            );
            assert!(
                f64::abs(eq4(d_3_0, f_3_0) - 0_f64) < eps,
                "result: {}, exp: {}",
                eq4(d_3_0, f_3_0),
                0_f64
            );
            //
            assert!(
                f64::abs(eq1(d_2_0, d_4_0, f_2_0, f_4_0) - 0_f64) < eps,
                "result: {}, exp: {}",
                eq1(d_2_0, d_4_0, f_2_0, f_4_0),
                0_f64
            );
            assert!(
                f64::abs(eq2(d_2_0, d_4_0, f_2_0, f_4_0) - (1_f64 + mu_0) * 4_f64 * alpha) < eps,
                "result: {}, exp: {}",
                eq2(d_2_0, d_4_0, f_2_0, f_4_0),
                (1_f64 + mu_0) * 4_f64 * alpha
            );
            assert!(
                f64::abs(eq3(d_2_0, d_4_0, f_2_0, f_4_0) - 0_f64) < eps,
                "result: {}, exp: {}",
                eq3(d_2_0, d_4_0, f_2_0, f_4_0),
                0_f64
            );
            assert!(
                f64::abs(eq4(d_4_0, f_4_0) - 0_f64) < eps,
                "result: {}, exp: {}",
                eq4(d_4_0, f_4_0),
                0_f64
            );
            //
            assert!(
                f64::abs(eq1(d_1_1, d_3_1, f_1_1, f_3_1) - 0_f64) < eps,
                "result: {}, exp: {}",
                eq1(d_1_1, d_3_1, f_1_1, f_3_1),
                0_f64
            );
            assert!(
                f64::abs(eq2(d_1_1, d_3_1, f_1_1, f_3_1) - 0_f64) < eps,
                "result: {}, exp: {}",
                eq2(d_1_1, d_3_1, f_1_1, f_3_1),
                0_f64
            );
            assert!(
                f64::abs(eq3(d_1_1, d_3_1, f_1_1, f_3_1) - (1_f64 + mu_0) * 4_f64 * alpha) < eps,
                "result: {}, exp: {}",
                eq3(d_1_1, d_3_1, f_1_1, f_3_1),
                (1_f64 + mu_0) * 4_f64 * alpha
            );
            assert!(
                f64::abs(eq4(d_3_1, f_3_1) - 0_f64) < eps,
                "result: {}, exp: {}",
                eq4(d_3_1, f_3_1),
                0_f64
            );
            //
            assert!(
                f64::abs(eq1(d_2_1, d_4_1, f_2_1, f_4_1) - 0_f64) < eps,
                "result: {}, exp: {}",
                eq1(d_2_1, d_4_1, f_2_1, f_4_1),
                0_f64
            );
            assert!(
                f64::abs(eq2(d_2_1, d_4_1, f_2_1, f_4_1) - 0_f64) < eps,
                "result: {}, exp: {}",
                eq2(d_2_1, d_4_1, f_2_1, f_4_1),
                0_f64
            );
            assert!(
                f64::abs(eq3(d_2_1, d_4_1, f_2_1, f_4_1) - 0_f64) < eps,
                "result: {}, exp: {}",
                eq3(d_2_1, d_4_1, f_2_1, f_4_1),
                0_f64
            );
            assert!(
                f64::abs(eq4(d_4_1, f_4_1) - (1_f64 + mu_0) * 4_f64 * alpha) < eps,
                "result: {}, exp: {}",
                eq4(d_4_1, f_4_1),
                (1_f64 + mu_0) * 4_f64 * alpha
            );
        }
    }
}
