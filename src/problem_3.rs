use crate::utils::{g, lambda, mu_0};
use std::f64::consts::PI;

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
