use crate::{
    integration::definite_integral,
    polynomials::chebyshev,
    utils::{g, lambda, mu_0, sum_calc},
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
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    m: usize,
    eps: f64,
) -> f64 {
    let coef = a * (f64::cosh(b * PI / a) - 1_f64);
    let f = |x| {
        let a1 = chebyshev(f64::cos(x * PI / (2_f64 * a)), 2 * m + 1);
        let a2 = f64::sqrt(1_f64 - x * x);
        let f = |i| {
            let alpha = PI / a * (i as f64 - 0.5);
            let pn = definite_integral(0_f64, a, 100, eps, &|x| {
                load_function(x) * f64::cos(alpha * x)
            });
            let ((_, psi_2_0, _, _), (_, _, _, _)) = psi_1(b, b, alpha, mu_0, g, lambda);
            let ((_, _, _, der_psi_4_0), (_, _, _, _)) = der_psi_1(b, b, alpha, mu_0, g, lambda);
            ((2_f64 * g + lambda) * der_psi_4_0 + lambda * psi_2_0 / alpha)
                * pn
                * f64::cos(alpha * x)
        };
        let a3 = sum_calc(0_f64, &f, eps, 1, 10) / a;
        let p = load_function(x);
        a1 * (a3 - p) / a2
    };
    definite_integral(-1_f64, 1_f64, 100, eps, &f) / coef
}

fn g_k<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    k: usize,
    eps: f64,
) -> f64 {
    let f = |x| {
        let a1 = chebyshev(f64::cosh(PI * b / a) * x + 1_f64, k);
        let a2 = f64::sqrt(1_f64 - x * x);

        let h = b - a * f64::acos(f64::cosh(PI * b / a) * x + 1_f64) / PI;
        let a3 = -a * g * (2_f64 * g + lambda) * b * mu_0 * mu_0 * h / (1_f64 + mu_0) / PI;
        let f1 = |i| {
            let alpha = PI / a * (i as f64 - 0.5);
            let ((_, _, _, _), (b_psi_1_1, b_psi_2_1, _, _)) = psi_1(x, b, alpha, mu_0, g, lambda);
            let ((_, _, _, _), (_, _, b_der_psi_3_1, b_der_psi_4_1)) =
                der_psi_1(x, b, alpha, mu_0, g, lambda);
            let ((h_psi_1_0, _, h_psi_3_0, _), (__, _, _, _)) = psi_1(h, b, alpha, mu_0, g, lambda);

            let a1 = f64::sin(alpha * (x + a)) + f64::sin(alpha * (a - x));
            let a2 = (2_f64 * g) * (b_der_psi_3_1 * h_psi_1_0 / alpha + b_der_psi_4_1 * h_psi_3_0)
                + lambda * (b_psi_1_1 * h_psi_1_0 / alpha / alpha + b_psi_2_1 * h_psi_3_0 / alpha);

            f64::exp(h - b) * a1 * a2 / alpha
        };
        let f2 = |i| {
            let i = i as f64;

            let a1 = f64::sin((2_f64 * i + 1_f64) * (PI / a / 2_f64) * (x + a))
                + f64::sin((2_f64 * i + 1_f64) * (PI / a / 2_f64) * (a - x));

            f64::exp((2_f64 * i + 1_f64) * (PI / a / 2_f64) * (h - b)) * a1 / (2_f64 * i + 1_f64)
        };
        let a4_1 = sum_calc(0_f64, &f1, eps, 1, 10);
        let a4_2 = sum_calc(0_f64, &f2, eps, 0, 10);
        let a4 = a4_1 - a3 * a4_2;

        2_f64 * a1 * a4 / a2 / a3
    };
    definite_integral(-1_f64, 1_f64, 100, eps, &f) / PI
}

fn psi_1(
    x: f64,
    b: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
) -> ((f64, f64, f64, f64), (f64, f64, f64, f64)) {
    let coef = (1_f64 + mu_0) * 4_f64;
    let e = f64::exp(-2_f64 * alpha * x);

    let (
        (d_1_0, d_2_0, d_3_0, d_4_0, f_1_0, f_2_0, f_3_0, f_4_0),
        (d_1_1, d_2_1, d_3_1, d_4_1, f_1_1, f_2_1, f_3_1, f_4_1),
    ) = coefficients_1(b, alpha, mu_0, g, lambda);

    let psi_1_0 = ((x * mu_0 + 2_f64 / alpha + mu_0 / alpha) * d_1_0 + x * mu_0 * d_3_0)
        + e * ((x * mu_0 - 2_f64 / alpha - mu_0 / alpha) * f_1_0 - x * mu_0 * f_3_0);
    let psi_2_0 = ((x * mu_0 + 2_f64 / alpha + mu_0 / alpha) * d_2_0 + x * mu_0 * d_4_0)
        + e * ((x * mu_0 - 2_f64 / alpha - mu_0 / alpha) * f_2_0 - x * mu_0 * f_4_0);
    let psi_3_0 = (-x * mu_0 * d_1_0 + (-x * mu_0 + 2_f64 / alpha + mu_0 / alpha) * d_3_0)
        + e * (x * mu_0 * f_1_0 + (-x * mu_0 - 2_f64 / alpha - mu_0 / alpha) * f_3_0);
    let psi_4_0 = (-x * mu_0 * d_2_0 + (-x * mu_0 + 2_f64 / alpha + mu_0 / alpha) * d_4_0)
        + e * (x * mu_0 * f_2_0 + (-x * mu_0 - 2_f64 / alpha - mu_0 / alpha) * f_4_0);

    let psi_1_1 = ((x * mu_0 + 2_f64 / alpha + mu_0 / alpha) * d_1_1 + x * mu_0 * d_3_1)
        + e * ((x * mu_0 - 2_f64 / alpha - mu_0 / alpha) * f_1_1 - x * mu_0 * f_3_1);
    let psi_2_1 = ((x * mu_0 + 2_f64 / alpha + mu_0 / alpha) * d_2_1 + x * mu_0 * d_4_1)
        + e * ((x * mu_0 - 2_f64 / alpha - mu_0 / alpha) * f_2_1 - x * mu_0 * f_4_1);
    let psi_3_1 = (-x * mu_0 * d_1_1 + (-x * mu_0 + 2_f64 / alpha + mu_0 / alpha) * d_3_1)
        + e * (x * mu_0 * f_1_1 + (-x * mu_0 - 2_f64 / alpha - mu_0 / alpha) * f_3_1);
    let psi_4_1 = (-x * mu_0 * d_2_1 + (-x * mu_0 + 2_f64 / alpha + mu_0 / alpha) * d_4_1)
        + e * (x * mu_0 * f_2_1 + (-x * mu_0 - 2_f64 / alpha - mu_0 / alpha) * f_4_1);

    (
        (
            psi_1_0 / coef,
            psi_2_0 / coef,
            psi_3_0 / coef,
            psi_4_0 / coef,
        ),
        (
            psi_1_1 / coef,
            psi_2_1 / coef,
            psi_3_1 / coef,
            psi_4_1 / coef,
        ),
    )
}

fn der_psi_1(
    x: f64,
    b: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
) -> ((f64, f64, f64, f64), (f64, f64, f64, f64)) {
    let coef = (1_f64 + mu_0) * 4_f64;
    let e = f64::exp(-2_f64 * alpha * x);

    let (
        (d_1_0, d_2_0, d_3_0, d_4_0, f_1_0, f_2_0, f_3_0, f_4_0),
        (d_1_1, d_2_1, d_3_1, d_4_1, f_1_1, f_2_1, f_3_1, f_4_1),
    ) = coefficients_1(b, alpha, mu_0, g, lambda);

    let psi_1_0 = ((x * mu_0 + 2_f64 / alpha + 2_f64 * mu_0 / alpha) * d_1_0
        + (x * mu_0 + mu_0 / alpha) * d_3_0)
        + e * ((-x * mu_0 + 2_f64 / alpha + 2_f64 * mu_0 / alpha) * f_1_0
            + (x * mu_0 - mu_0 / alpha) * f_3_0);
    let psi_2_0 = ((x * mu_0 + 2_f64 / alpha + 2_f64 * mu_0 / alpha) * d_2_0
        + (x * mu_0 + mu_0 / alpha) * d_4_0)
        + e * ((-x * mu_0 + 2_f64 / alpha + 2_f64 * mu_0 / alpha) * f_2_0
            + (x * mu_0 - mu_0 / alpha) * f_4_0);
    let psi_3_0 = ((-x * mu_0 - mu_0 / alpha) * d_1_0 + (-x * mu_0 + 2_f64 / alpha) * d_3_0)
        + e * ((-x * mu_0 + mu_0 / alpha) * f_1_0 + (x * mu_0 + 2_f64 / alpha) * f_3_0);
    let psi_4_0 = ((-x * mu_0 - mu_0 / alpha) * d_2_0 + (-x * mu_0 + 2_f64 / alpha) * d_4_0)
        + e * ((-x * mu_0 + mu_0 / alpha) * f_2_0 + (x * mu_0 + 2_f64 / alpha) * f_4_0);

    let psi_1_1 = ((x * mu_0 + 2_f64 / alpha + 2_f64 * mu_0 / alpha) * d_1_1
        + (x * mu_0 + mu_0 / alpha) * d_3_1)
        + e * ((-x * mu_0 + 2_f64 / alpha + 2_f64 * mu_0 / alpha) * f_1_1
            + (x * mu_0 - mu_0 / alpha) * f_3_1);
    let psi_2_1 = ((x * mu_0 + 2_f64 / alpha + 2_f64 * mu_0 / alpha) * d_2_1
        + (x * mu_0 + mu_0 / alpha) * d_4_1)
        + e * ((-x * mu_0 + 2_f64 / alpha + 2_f64 * mu_0 / alpha) * f_2_1
            + (x * mu_0 - mu_0 / alpha) * f_4_1);
    let psi_3_1 = ((-x * mu_0 - mu_0 / alpha) * d_1_1 + (-x * mu_0 + 2_f64 / alpha) * d_3_1)
        + e * ((-x * mu_0 + mu_0 / alpha) * f_1_1 + (x * mu_0 + 2_f64 / alpha) * f_3_1);
    let psi_4_1 = ((-x * mu_0 - mu_0 / alpha) * d_2_1 + (-x * mu_0 + 2_f64 / alpha) * d_4_1)
        + e * ((-x * mu_0 + mu_0 / alpha) * f_2_1 + (x * mu_0 + 2_f64 / alpha) * f_4_1);

    (
        (
            psi_1_0 / coef,
            psi_2_0 / coef,
            psi_3_0 / coef,
            psi_4_0 / coef,
        ),
        (
            psi_1_1 / coef,
            psi_2_1 / coef,
            psi_3_1 / coef,
            psi_4_1 / coef,
        ),
    )
}

fn psi_2(
    x: f64,
    b: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
) -> ((f64, f64, f64, f64), (f64, f64, f64, f64)) {
    let e1 = f64::exp(alpha * x);
    let e2 = f64::exp(-alpha * b);

    let ((psi_1_0, psi_2_0, psi_3_0, psi_4_0), (psi_1_1, psi_2_1, psi_3_1, psi_4_1)) =
        psi_1(x, b, alpha, mu_0, g, lambda);

    let psi_1_0 = e1 * e2 * psi_1_0 / alpha;
    let psi_2_0 = e1 * e2 * psi_2_0 / alpha;
    let psi_3_0 = e1 * e2 * psi_3_0 / alpha;
    let psi_4_0 = e1 * e2 * psi_4_0 / alpha;

    let psi_1_1 = e1 * psi_1_1;
    let psi_2_1 = alpha * e1 * psi_2_1;
    let psi_3_1 = e1 * psi_3_1;
    let psi_4_1 = alpha * e1 * psi_4_1;

    (
        (psi_1_0, psi_2_0, psi_3_0, psi_4_0),
        (psi_1_1, psi_2_1, psi_3_1, psi_4_1),
    )
}

fn der_psi_2(
    x: f64,
    b: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
) -> ((f64, f64, f64, f64), (f64, f64, f64, f64)) {
    let e1 = f64::exp(alpha * x);
    let e2 = f64::exp(-alpha * b);

    let ((psi_1_0, psi_2_0, psi_3_0, psi_4_0), (psi_1_1, psi_2_1, psi_3_1, psi_4_1)) =
        der_psi_1(x, b, alpha, mu_0, g, lambda);

    let psi_1_0 = e1 * e2 * psi_1_0;
    let psi_2_0 = e1 * e2 * psi_2_0;
    let psi_3_0 = e1 * e2 * psi_3_0;
    let psi_4_0 = e1 * e2 * psi_4_0;

    let psi_1_1 = alpha * e1 * psi_1_1;
    let psi_2_1 = alpha * alpha * e1 * psi_2_1;
    let psi_3_1 = alpha * e1 * psi_3_1;
    let psi_4_1 = alpha * alpha * e1 * psi_4_1;

    (
        (psi_1_0, psi_2_0, psi_3_0, psi_4_0),
        (psi_1_1, psi_2_1, psi_3_1, psi_4_1),
    )
}

fn coefficients_1(
    b: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
) -> (
    (f64, f64, f64, f64, f64, f64, f64, f64),
    (f64, f64, f64, f64, f64, f64, f64, f64),
) {
    let coef = (1_f64 + mu_0) * 4_f64;
    let e = f64::exp(-alpha * b);

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

    let y1 = x2 + e * e * x4;
    let y2 = e * e * x3 - x1;
    let y3 = x6 + e * e * x8;
    let y4 = e * e * x7 - x5;

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
        * ((y2 * (x1 * y4 - x5 * y2) * (x10 * x12 + x10 * x11)
            + x11 * x9 * (x1 * y4 - x5 * y2) * y1
            + x11 * x9 * x1 * (y3 * y2 - y1 * y4))
            / (x11 * x9 * x9 * y2 * (y3 * y2 - y1 * y4)))
        + coef / x9;
    let d_3_1 = -coef * (x12 * (x1 * y4 - x5 * y2) / (x11 * x9 * (y3 * y2 - y1 * y4)));
    let f_1_1 = -coef
        * ((x1 * (y3 * y2 - y1 * y4) + (x1 * y4 - x5 * y2) * y1) / (x9 * y2 * (y3 * y2 - y1 * y4)));
    let f_3_1 = coef * (x1 * y4 - x5 * y2) / (x9 * (y3 * y2 - y1 * y4));

    let d_2_1 = coef
        * (((x5 * x10 - x1 * x9) * y2 - (x1 * x10 - x2 * x9) * y4) * (x10 * x11 + x10 * x12)
            / (x9 * x9 * x11 * x11 * (y2 * y3 - y1 * y4))
            - ((x1 * x10 - x2 * x9) * y3 - (x5 * x10 - x6 * x9) * y1)
                / (x9 * x11 * (y2 * y3 - y1 * y4)))
        - coef * x10 / (x9 * x11);
    let d_4_1 = -coef * ((x5 * x10 - x6 * x9) * y2 - (x1 * x10 - x2 * x9) * y4) * x12
        / (x9 * x11 * x11 * (y2 * y3 - y1 * y4))
        + coef / x11;
    let f_2_1 = coef * ((x1 * x10 - x2 * x9) * y3 - (x5 * x10 - x6 * x9) * y1)
        / (x9 * x11 * (y2 * y3 - y1 * y4));
    let f_4_1 = coef * ((x5 * x10 - x6 * x9) * y2 - (x1 * x10 - x2 * x9) * y4)
        / (x9 * x11 * (y2 * y3 - y1 * y4));

    (
        (d_1_0, d_2_0, d_3_0, d_4_0, f_1_0, f_2_0, f_3_0, f_4_0),
        (d_1_1, d_2_1, d_3_1, d_4_1, f_1_1, f_2_1, f_3_1, f_4_1),
    )
}

fn coefficients_2(
    b: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
) -> (
    (f64, f64, f64, f64, f64, f64, f64, f64),
    (f64, f64, f64, f64, f64, f64, f64, f64),
) {
    let e = f64::exp(-alpha * b);

    let (
        (d_1_0, d_2_0, d_3_0, d_4_0, f_1_0, f_2_0, f_3_0, f_4_0),
        (d_1_1, d_2_1, d_3_1, d_4_1, f_1_1, f_2_1, f_3_1, f_4_1),
    ) = coefficients_1(b, alpha, mu_0, g, lambda);

    let d_1_0 = e * d_1_0 / alpha;
    let d_2_0 = e * d_2_0 / alpha;
    let d_3_0 = e * d_3_0 / alpha;
    let d_4_0 = e * d_4_0 / alpha;

    let f_1_0 = e * f_1_0 / alpha;
    let f_2_0 = e * f_2_0 / alpha;
    let f_3_0 = e * f_3_0 / alpha;
    let f_4_0 = e * f_4_0 / alpha;

    let d_1_1 = d_1_1;
    let d_2_1 = d_2_1 * alpha;
    let d_3_1 = d_3_1;
    let d_4_1 = d_4_1 * alpha;

    let f_1_1 = f_1_1;
    let f_2_1 = f_2_1 * alpha;
    let f_3_1 = f_3_1;
    let f_4_1 = f_4_1 * alpha;

    (
        (d_1_0, d_2_0, d_3_0, d_4_0, f_1_0, f_2_0, f_3_0, f_4_0),
        (d_1_1, d_2_1, d_3_1, d_4_1, f_1_1, f_2_1, f_3_1, f_4_1),
    )
}

fn unknown_function(x: f64, a: f64, b: f64, mu_0: f64, g: f64, lambda: f64, eps: f64) -> f64 {
    let h = b - a * f64::acos(f64::cosh(PI * b / a) * x + 1_f64) / PI;
    let a1 = -a * g * (2_f64 * g + lambda) * b * mu_0 * mu_0 * h / (1_f64 + mu_0) / PI;
    let a2 = 2_f64
        * f64::sqrt(
            (f64::cosh(PI * b / a) * h + 1_f64) * (f64::cosh(PI * b / a) * h + 1_f64) - 1_f64,
        );
    let a3 = f64::sqrt(1_f64 - h * h);
    let f = |i| chebyshev(f64::cosh(PI * b / a * h + 1_f64), 2 * i + 1);
    let a4 = sum_calc(0_f64, &f, eps, 0, 10);

    a2 * a4 / a1 / a3
}

#[cfg(test)]
mod tests {
    use super::*;

    fn der_2_psi_1(
        x: f64,
        b: f64,
        alpha: f64,
        mu_0: f64,
        g: f64,
        lambda: f64,
    ) -> ((f64, f64, f64, f64), (f64, f64, f64, f64)) {
        let coef = (1_f64 + mu_0) * 4_f64;
        let e = f64::exp(-2_f64 * alpha * x);

        let (
            (d_1_0, d_2_0, d_3_0, d_4_0, f_1_0, f_2_0, f_3_0, f_4_0),
            (d_1_1, d_2_1, d_3_1, d_4_1, f_1_1, f_2_1, f_3_1, f_4_1),
        ) = coefficients_1(b, alpha, mu_0, g, lambda);

        let psi_1_0 = ((x * mu_0 + 2_f64 / alpha + 3_f64 * mu_0 / alpha) * d_1_0
            + (x * mu_0 + 2_f64 * mu_0 / alpha) * d_3_0)
            + e * ((x * mu_0 - 2_f64 / alpha - 3_f64 * mu_0 / alpha) * f_1_0
                + (-x * mu_0 + 2_f64 * mu_0 / alpha) * f_3_0);
        let psi_2_0 = ((x * mu_0 + 2_f64 / alpha + 3_f64 * mu_0 / alpha) * d_2_0
            + (x * mu_0 + 2_f64 * mu_0 / alpha) * d_4_0)
            + e * ((x * mu_0 - 2_f64 / alpha - 3_f64 * mu_0 / alpha) * f_2_0
                + (-x * mu_0 + 2_f64 * mu_0 / alpha) * f_4_0);
        let psi_3_0 = ((-x * mu_0 - 2_f64 * mu_0 / alpha) * d_1_0
            + (-x * mu_0 + 2_f64 / alpha - mu_0 / alpha) * d_3_0)
            + e * ((x * mu_0 - 2_f64 * mu_0 / alpha) * f_1_0
                + (-x * mu_0 - 2_f64 / alpha + mu_0 / alpha) * f_3_0);
        let psi_4_0 = ((-x * mu_0 - 2_f64 * mu_0 / alpha) * d_2_0
            + (-x * mu_0 + 2_f64 / alpha - mu_0 / alpha) * d_4_0)
            + e * ((x * mu_0 - 2_f64 * mu_0 / alpha) * f_2_0
                + (-x * mu_0 - 2_f64 / alpha + mu_0 / alpha) * f_4_0);

        let psi_1_1 = ((x * mu_0 + 2_f64 / alpha + 3_f64 * mu_0 / alpha) * d_1_1
            + (x * mu_0 + 2_f64 * mu_0 / alpha) * d_3_1)
            + e * ((x * mu_0 - 2_f64 / alpha - 3_f64 * mu_0 / alpha) * f_1_1
                + (-x * mu_0 + 2_f64 * mu_0 / alpha) * f_3_1);
        let psi_2_1 = ((x * mu_0 + 2_f64 / alpha + 3_f64 * mu_0 / alpha) * d_2_1
            + (x * mu_0 + 2_f64 * mu_0 / alpha) * d_4_1)
            + e * ((x * mu_0 - 2_f64 / alpha - 3_f64 * mu_0 / alpha) * f_2_1
                + (-x * mu_0 + 2_f64 * mu_0 / alpha) * f_4_1);
        let psi_3_1 = ((-x * mu_0 - 2_f64 * mu_0 / alpha) * d_1_1
            + (-x * mu_0 + 2_f64 / alpha - mu_0 / alpha) * d_3_1)
            + e * ((x * mu_0 - 2_f64 * mu_0 / alpha) * f_1_1
                + (-x * mu_0 - 2_f64 / alpha + mu_0 / alpha) * f_3_1);
        let psi_4_1 = ((-x * mu_0 - 2_f64 * mu_0 / alpha) * d_2_1
            + (-x * mu_0 + 2_f64 / alpha - mu_0 / alpha) * d_4_1)
            + e * ((x * mu_0 - 2_f64 * mu_0 / alpha) * f_2_1
                + (-x * mu_0 - 2_f64 / alpha + mu_0 / alpha) * f_4_1);

        (
            (
                psi_1_0 / coef,
                psi_2_0 / coef,
                psi_3_0 / coef,
                psi_4_0 / coef,
            ),
            (
                psi_1_1 / coef,
                psi_2_1 / coef,
                psi_3_1 / coef,
                psi_4_1 / coef,
            ),
        )
    }

    fn der_2_psi_2(
        x: f64,
        b: f64,
        alpha: f64,
        mu_0: f64,
        g: f64,
        lambda: f64,
    ) -> ((f64, f64, f64, f64), (f64, f64, f64, f64)) {
        let e1 = f64::exp(alpha * x);
        let e2 = f64::exp(-alpha * b);

        let ((psi_1_0, psi_2_0, psi_3_0, psi_4_0), (psi_1_1, psi_2_1, psi_3_1, psi_4_1)) =
            der_2_psi_1(x, b, alpha, mu_0, g, lambda);

        let psi_1_0 = alpha * e1 * e2 * psi_1_0;
        let psi_2_0 = alpha * e1 * e2 * psi_2_0;
        let psi_3_0 = alpha * e1 * e2 * psi_3_0;
        let psi_4_0 = alpha * e1 * e2 * psi_4_0;

        let psi_1_1 = alpha * alpha * e1 * psi_1_1;
        let psi_2_1 = alpha * alpha * alpha * e1 * psi_2_1;
        let psi_3_1 = alpha * alpha * e1 * psi_3_1;
        let psi_4_1 = alpha * alpha * alpha * e1 * psi_4_1;

        (
            (psi_1_0, psi_2_0, psi_3_0, psi_4_0),
            (psi_1_1, psi_2_1, psi_3_1, psi_4_1),
        )
    }

    #[test]
    fn coefficients_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;
        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);

        for i in 1..5 {
            println!("-----");
            let alpha = PI / a * (i as f64 - 0.5);

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
            ) = coefficients_2(b, alpha, mu_0, g, lambda);

            let eq1 = |d1, d2, f1, f2| e1 * x1 * d1 + e1 * x2 * d2 + e2 * x3 * f1 + e2 * x4 * f2;
            let eq2 = |d1, d2, f1, f2| e1 * x5 * d1 + e1 * x6 * d2 + e2 * x7 * f1 + e2 * x8 * f2;
            let eq3 = |d1, d2, f1, f2| x9 * d1 + x10 * d2 + x9 * f1 - x10 * f2;
            let eq4 = |d1, f1| x11 * d1 + x12 * f1;
            //
            println!(
                "left: {0}, right: {1}",
                eq1(d_1_0, d_3_0, f_1_0, f_3_0),
                (1_f64 + mu_0) * 4_f64 * alpha
            );
            println!(
                "left: {0}, right: {1}",
                eq2(d_1_0, d_3_0, f_1_0, f_3_0),
                0_f64
            );
            println!(
                "left: {0}, right: {1}",
                eq2(d_1_0, d_3_0, f_1_0, f_3_0),
                0_f64
            );
            println!("left: {0}, right: {1}", eq4(d_3_0, f_3_0), 0_f64);
            //
            println!(
                "left: {0}, right: {1}",
                eq1(d_2_0, d_4_0, f_2_0, f_4_0),
                0_f64
            );
            println!(
                "left: {0}, right: {1}",
                eq2(d_2_0, d_4_0, f_2_0, f_4_0),
                (1_f64 + mu_0) * 4_f64 * alpha
            );
            println!(
                "left: {0}, right: {1}",
                eq3(d_2_0, d_4_0, f_2_0, f_4_0),
                0_f64
            );
            println!("left: {0}, right: {1}", eq4(d_4_0, f_4_0), 0_f64);
            //
            println!(
                "left: {0}, right: {1}",
                eq1(d_1_1, d_3_1, f_1_1, f_3_1),
                0_f64
            );
            println!(
                "left: {0}, right: {1}",
                eq2(d_1_1, d_3_1, f_1_1, f_3_1),
                0_f64
            );
            println!(
                "left: {0}, right: {1}",
                eq3(d_1_1, d_3_1, f_1_1, f_3_1),
                (1_f64 + mu_0) * 4_f64 * alpha
            );
            println!("left: {0}, right: {1}", eq4(d_3_1, f_3_1), 0_f64);
            //
            println!(
                "left: {0}, right: {1}",
                eq1(d_2_1, d_4_1, f_2_1, f_4_1),
                0_f64
            );
            println!(
                "left: {0}, right: {1}",
                eq2(d_2_1, d_4_1, f_2_1, f_4_1),
                0_f64
            );
            println!(
                "left: {0}, right: {1}",
                eq3(d_2_1, d_4_1, f_2_1, f_4_1),
                0_f64
            );
            println!(
                "left: {0}, right: {1}",
                eq4(d_4_1, f_4_1),
                (1_f64 + mu_0) * 4_f64 * alpha
            );
        }
    }

    #[test]
    fn equation_psi_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;

        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);

        for x in 0..(b - 1_f64) as usize {
            let x = x as f64;
            for i in 1..20 {
                println!("----");
                let alpha = PI / a * (i as f64 - 0.5);

                let ((psi_1_0, psi_2_0, psi_3_0, psi_4_0), (psi_1_1, psi_2_1, psi_3_1, psi_4_1)) =
                    psi_2(x, b, alpha, mu_0, g, lambda);
                let (
                    (der_psi_1_0, der_psi_2_0, der_psi_3_0, der_psi_4_0),
                    (der_psi_1_1, der_psi_2_1, der_psi_3_1, der_psi_4_1),
                ) = der_psi_2(x, b, alpha, mu_0, g, lambda);
                let (
                    (der_2_psi_1_0, der_2_psi_2_0, der_2_psi_3_0, der_2_psi_4_0),
                    (der_2_psi_1_1, der_2_psi_2_1, der_2_psi_3_1, der_2_psi_4_1),
                ) = der_2_psi_2(x, b, alpha, mu_0, g, lambda);

                let line_1 = der_2_psi_1_0
                    + (-alpha * mu_0) * der_psi_3_0
                    + (-alpha * alpha - alpha * alpha * mu_0) * psi_1_0;
                let line_2 = der_2_psi_2_0
                    + (-alpha * mu_0) * der_psi_4_0
                    + (-alpha * alpha - alpha * alpha * mu_0) * psi_2_0;
                let line_3 = der_2_psi_3_0 * (1_f64 + mu_0)
                    + (alpha * mu_0) * der_psi_1_0
                    + (-alpha * alpha) * psi_3_0;
                let line_4 = der_2_psi_4_0 * (1_f64 + mu_0)
                    + (alpha * mu_0) * der_psi_2_0
                    + (-alpha * alpha) * psi_4_0;

                println!("0 line 1: {line_1}");
                println!("0 line 2: {line_2}");
                println!("0 line 3: {line_3}");
                println!("0 line 4: {line_4}");

                let line_1 = der_2_psi_1_1
                    + (-alpha * mu_0) * der_psi_3_1
                    + (-alpha * alpha - alpha * alpha * mu_0) * psi_1_1;
                let line_2 = der_2_psi_2_1
                    + (-alpha * mu_0) * der_psi_4_1
                    + (-alpha * alpha - alpha * alpha * mu_0) * psi_2_1;
                let line_3 = der_2_psi_3_1 * (1_f64 + mu_0)
                    + (alpha * mu_0) * der_psi_1_1
                    + (-alpha * alpha) * psi_3_1;
                let line_4 = der_2_psi_4_1 * (1_f64 + mu_0)
                    + (alpha * mu_0) * der_psi_2_1
                    + (-alpha * alpha) * psi_4_1;

                println!("1 line 1: {line_1}");
                println!("1 line 2: {line_2}");
                println!("1 line 3: {line_3}");
                println!("1 line 4: {line_4}");
            }
        }
    }

    #[test]
    fn boundary_condition_1_psi_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;

        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);

        for i in 1..20 {
            println!("----");
            let alpha = PI / a * (i as f64 - 0.5);

            let ((psi_1_0, psi_2_0, psi_3_0, psi_4_0), (psi_1_1, psi_2_1, psi_3_1, psi_4_1)) =
                psi_2(b, b, alpha, mu_0, g, lambda);
            let (
                (der_psi_1_0, der_psi_2_0, der_psi_3_0, der_psi_4_0),
                (der_psi_1_1, der_psi_2_1, der_psi_3_1, der_psi_4_1),
            ) = der_psi_2(b, b, alpha, mu_0, g, lambda);

            let line_1 = der_psi_1_0 + (-alpha) * psi_3_0;
            let line_2 = der_psi_2_0 + (-alpha) * psi_4_0;
            let line_3 = (2_f64 * g + lambda) * der_psi_3_0 + (alpha * lambda) * psi_1_0;
            let line_4 = (2_f64 * g + lambda) * der_psi_4_0 + (alpha * lambda) * psi_2_0;

            println!("0 line 1: {line_1}");
            println!("0 line 2: {line_2}");
            println!("0 line 3: {line_3}");
            println!("0 line 4: {line_4}");

            let line_1 = der_psi_1_1 + (-alpha) * psi_3_1;
            let line_2 = der_psi_2_1 + (-alpha) * psi_4_1;
            let line_3 = (2_f64 * g + lambda) * der_psi_3_1 + (alpha * lambda) * psi_1_1;
            let line_4 = (2_f64 * g + lambda) * der_psi_4_1 + (alpha * lambda) * psi_2_1;

            println!("1 line 1: {line_1}");
            println!("1 line 2: {line_2}");
            println!("1 line 3: {line_3}");
            println!("1 line 4: {line_4}");
        }
    }

    #[test]
    fn boundary_condition_2_psi_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;

        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);

        for i in 1..20 {
            println!("----");
            let alpha = PI / a * (i as f64 - 0.5);

            let ((_, _, psi_3_0, psi_4_0), (_, _, psi_3_1, psi_4_1)) =
                psi_2(0_f64, b, alpha, mu_0, g, lambda);
            let ((der_psi_1_0, der_psi_2_0, _, _), (der_psi_1_1, der_psi_2_1, _, _)) =
                der_psi_2(0_f64, b, alpha, mu_0, g, lambda);

            let line_1 = der_psi_1_0 + (-alpha) * psi_3_0;
            let line_2 = der_psi_2_0 + (-alpha) * psi_4_0;
            let line_3 = psi_3_0;
            let line_4 = psi_4_0;

            println!("0 line 1: {line_1}");
            println!("0 line 2: {line_2}");
            println!("0 line 3: {line_3}");
            println!("0 line 4: {line_4}");

            let line_1 = der_psi_1_1 + (-alpha) * psi_3_1;
            let line_2 = der_psi_2_1 + (-alpha) * psi_4_1;
            let line_3 = psi_3_1;
            let line_4 = psi_4_1;

            println!("1 line 1: {line_1}");
            println!("1 line 2: {line_2}");
            println!("1 line 3: {line_3}");
            println!("1 line 4: {line_4}");
        }
    }
}
