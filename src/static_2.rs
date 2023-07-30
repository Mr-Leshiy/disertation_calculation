use crate::{
    integration::{
        definite_integral, definite_integral_limit, sqrt_gauss_integral, sqrt_gauss_integral_finit,
    },
    polynomials::chebyshev,
    utils::sum_calc_finit,
};
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use std::f64::consts::PI;

type Matrix<const N: usize, const M: usize> = nalgebra::Matrix<
    f64,
    nalgebra::Const<N>,
    nalgebra::Const<M>,
    nalgebra::ArrayStorage<f64, N, M>,
>;

fn alpha_calc(a: f64, n: f64) -> f64 {
    PI / a * (n - 0.5)
}

fn a_functions(
    alpha: f64,
    mu_0: f64,
) -> (
    impl Fn(f64) -> f64,
    impl Fn(f64) -> f64,
    impl Fn(f64) -> f64,
    impl Fn(f64) -> f64,
    impl Fn(f64) -> f64,
    impl Fn(f64) -> f64,
    impl Fn(f64) -> f64,
    impl Fn(f64) -> f64,
    impl Fn(f64) -> f64,
    impl Fn(f64) -> f64,
    impl Fn(f64) -> f64,
    impl Fn(f64) -> f64,
    impl Fn(f64) -> f64,
    impl Fn(f64) -> f64,
    impl Fn(f64) -> f64,
    impl Fn(f64) -> f64,
) {
    let coef = (1.0 + mu_0) * 4.0 * alpha;
    let a_1 = move |y: f64| f64::exp(alpha * y) * (y * alpha * mu_0 + 2.0 + mu_0) / coef;
    let a_2 = move |y: f64| f64::exp(alpha * y) * (y * alpha * mu_0) / coef;
    let a_3 = move |y: f64| f64::exp(-alpha * y) * (y * alpha * mu_0 - 2.0 - mu_0) / coef;
    let a_4 = move |y: f64| f64::exp(-alpha * y) * (-y * alpha * mu_0) / coef;
    let a_5 = move |y: f64| f64::exp(alpha * y) * (-y * alpha * mu_0) / coef;
    let a_6 = move |y: f64| f64::exp(alpha * y) * (-y * alpha * mu_0 + 2.0 + mu_0) / coef;
    let a_7 = move |y: f64| f64::exp(-alpha * y) * (y * alpha * mu_0) / coef;
    let a_8 = move |y: f64| f64::exp(-alpha * y) * (-y * alpha * mu_0 - 2.0 - mu_0) / coef;

    let a_9 =
        move |y: f64| f64::exp(alpha * y) * (y * alpha * mu_0 + 2.0 + 2.0 * mu_0) * alpha / coef;
    let a_10 = move |y: f64| f64::exp(alpha * y) * (y * alpha * mu_0 + mu_0) * alpha / coef;
    let a_11 =
        move |y: f64| f64::exp(-alpha * y) * (-y * alpha * mu_0 + 2.0 + 2.0 * mu_0) * alpha / coef;
    let a_12 = move |y: f64| f64::exp(-alpha * y) * (y * alpha * mu_0 - mu_0) * alpha / coef;
    let a_13 = move |y: f64| f64::exp(alpha * y) * (-y * alpha * mu_0 - mu_0) * alpha / coef;
    let a_14 = move |y: f64| f64::exp(alpha * y) * (-y * alpha * mu_0 + 2.0) * alpha / coef;
    let a_15 = move |y: f64| f64::exp(-alpha * y) * (-y * alpha * mu_0 + mu_0) * alpha / coef;
    let a_16 = move |y: f64| f64::exp(-alpha * y) * (y * alpha * mu_0 + 2.0) * alpha / coef;

    (
        a_1, a_2, a_3, a_4, a_5, a_6, a_7, a_8, a_9, a_10, a_11, a_12, a_13, a_14, a_15, a_16,
    )
}

fn b_values(
    b: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
) -> (
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
    f64,
) {
    let (a_1, a_2, a_3, a_4, a_5, a_6, a_7, a_8, a_9, a_10, a_11, a_12, a_13, a_14, a_15, a_16) =
        a_functions(alpha, mu_0);

    let b_1 = a_9(b) - alpha * a_5(b);
    let b_2 = a_10(b) - alpha * a_6(b);
    let b_3 = a_11(b) - alpha * a_7(b);
    let b_4 = a_12(b) - alpha * a_8(b);
    let b_5 = (2.0 * g + lambda) * a_13(b) + alpha * lambda * a_1(b);
    let b_6 = (2.0 * g + lambda) * a_14(b) + alpha * lambda * a_2(b);
    let b_7 = (2.0 * g + lambda) * a_15(b) + alpha * lambda * a_3(b);
    let b_8 = (2.0 * g + lambda) * a_16(b) + alpha * lambda * a_4(b);
    let b_9 = a_9(0.0) - alpha * a_5(0.0);
    let b_10 = a_10(0.0) - alpha * a_6(0.0);
    let b_11 = a_11(0.0) - alpha * a_7(0.0);
    let b_12 = a_12(0.0) - alpha * a_8(0.0);
    let b_13 = a_5(0.0);
    let b_14 = a_6(0.0);
    let b_15 = a_7(0.0);
    let b_16 = a_8(0.0);

    (
        b_1, b_2, b_3, b_4, b_5, b_6, b_7, b_8, b_9, b_10, b_11, b_12, b_13, b_14, b_15, b_16,
    )
}

fn coefficients(
    b: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
) -> (
    ((f64, f64, f64, f64), (f64, f64, f64, f64)),
    ((f64, f64, f64, f64), (f64, f64, f64, f64)),
) {
    let (b_1, b_2, b_3, b_4, b_5, b_6, b_7, b_8, b_9, b_10, b_11, b_12, _b_13, b_14, _b_15, b_16) =
        b_values(b, alpha, mu_0, g, lambda);

    let dem = (b_8 * b_14 - b_6 * b_16) * (b_3 * b_9 - b_1 * b_11)
        + (b_10 * b_16 - b_12 * b_14) * (b_3 * b_5 - b_1 * b_7)
        + (b_2 * b_16 - b_4 * b_14) * (b_7 * b_9 - b_5 * b_11);

    let d_1_0 = (b_16 * (b_6 * b_11 - b_7 * b_10) + b_14 * (b_7 * b_12 - b_8 * b_11)) / dem;
    let d_2_0 = (b_16 * (b_3 * b_10 - b_2 * b_11) + b_14 * (b_4 * b_11 - b_3 * b_12)) / dem;
    let d_3_0 = b_16 * (b_7 * b_9 - b_5 * b_11) / dem;
    let d_4_0 = b_16 * (b_1 * b_11 - b_3 * b_9) / dem;

    let f_1_0 = (b_16 * (b_5 * b_10 - b_6 * b_9) + b_14 * (b_8 * b_9 - b_5 * b_12)) / dem;
    let f_2_0 = (b_16 * (b_2 * b_9 - b_1 * b_10) + b_14 * (b_1 * b_12 - b_4 * b_9)) / dem;
    let f_3_0 = b_14 * (b_11 * b_5 - b_7 * b_9) / dem;
    let f_4_0 = b_14 * (b_3 * b_9 - b_1 * b_11) / dem;

    let d_1_1 = (b_14 * b_3 * b_8 - b_14 * b_4 * b_7 + b_16 * b_2 * b_7 - b_16 * b_3 * b_6) / dem;
    let d_2_1 = (-b_10 * b_3 * b_8 + b_10 * b_4 * b_7 + b_11 * b_2 * b_8
        - b_11 * b_4 * b_6
        - b_12 * b_2 * b_7
        + b_12 * b_3 * b_6)
        / dem;
    let d_3_1 = (-b_1 * b_16 * b_7 + b_16 * b_3 * b_5) / dem;
    let d_4_1 = (-b_1 * b_11 * b_8 + b_1 * b_12 * b_7 + b_11 * b_4 * b_5 - b_12 * b_3 * b_5
        + b_3 * b_8 * b_9
        - b_4 * b_7 * b_9)
        / dem;

    let f_1_1 = (-b_1 * b_14 * b_8 + b_1 * b_16 * b_6 + b_14 * b_4 * b_5 - b_16 * b_2 * b_5) / dem;
    let f_2_1 = (b_1 * b_10 * b_8 - b_1 * b_12 * b_6 - b_10 * b_4 * b_5 + b_12 * b_2 * b_5
        - b_2 * b_8 * b_9
        + b_4 * b_6 * b_9)
        / dem;
    let f_3_1 = (b_1 * b_14 * b_7 - b_14 * b_3 * b_5) / dem;
    let f_4_1 = (-b_1 * b_10 * b_7 + b_1 * b_11 * b_6 + b_10 * b_3 * b_5 - b_11 * b_2 * b_5
        + b_2 * b_7 * b_9
        - b_3 * b_6 * b_9)
        / dem;

    (
        ((d_1_0, d_2_0, d_3_0, d_4_0), (f_1_0, f_2_0, f_3_0, f_4_0)),
        ((d_1_1, d_2_1, d_3_1, d_4_1), (f_1_1, f_2_1, f_3_1, f_4_1)),
    )
}

fn psi(
    y: f64,
    b: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
) -> ((f64, f64, f64, f64), (f64, f64, f64, f64)) {
    let (a_1, a_2, a_3, a_4, a_5, a_6, a_7, a_8, _, _, _, _, _, _, _, _) = a_functions(alpha, mu_0);
    let (
        ((d_1_0, d_2_0, d_3_0, d_4_0), (f_1_0, f_2_0, f_3_0, f_4_0)),
        ((d_1_1, d_2_1, d_3_1, d_4_1), (f_1_1, f_2_1, f_3_1, f_4_1)),
    ) = coefficients(b, alpha, mu_0, g, lambda);

    let psi_1_0 = a_1(y) * d_1_0 + a_2(y) * d_3_0 + a_3(y) * f_1_0 + a_4(y) * f_3_0;
    let psi_2_0 = a_1(y) * d_2_0 + a_2(y) * d_4_0 + a_3(y) * f_2_0 + a_4(y) * f_4_0;
    let psi_3_0 = a_5(y) * d_1_0 + a_6(y) * d_3_0 + a_7(y) * f_1_0 + a_8(y) * f_3_0;
    let psi_4_0 = a_5(y) * d_2_0 + a_6(y) * d_4_0 + a_7(y) * f_2_0 + a_8(y) * f_4_0;

    let psi_1_1 = a_1(y) * d_1_1 + a_2(y) * d_3_1 + a_3(y) * f_1_1 + a_4(y) * f_3_1;
    let psi_2_1 = a_1(y) * d_2_1 + a_2(y) * d_4_1 + a_3(y) * f_2_1 + a_4(y) * f_4_1;
    let psi_3_1 = a_5(y) * d_1_1 + a_6(y) * d_3_1 + a_7(y) * f_1_1 + a_8(y) * f_3_1;
    let psi_4_1 = a_5(y) * d_2_1 + a_6(y) * d_4_1 + a_7(y) * f_2_1 + a_8(y) * f_4_1;

    (
        (psi_1_0, psi_2_0, psi_3_0, psi_4_0),
        (psi_1_1, psi_2_1, psi_3_1, psi_4_1),
    )
}

fn der_psi(
    y: f64,
    b: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
) -> ((f64, f64, f64, f64), (f64, f64, f64, f64)) {
    let (_, _, _, _, _, _, _, _, a_9, a_10, a_11, a_12, a_13, a_14, a_15, a_16) =
        a_functions(alpha, mu_0);
    let (
        ((d_1_0, d_2_0, d_3_0, d_4_0), (f_1_0, f_2_0, f_3_0, f_4_0)),
        ((d_1_1, d_2_1, d_3_1, d_4_1), (f_1_1, f_2_1, f_3_1, f_4_1)),
    ) = coefficients(b, alpha, mu_0, g, lambda);

    let der_psi_1_0 = a_9(y) * d_1_0 + a_10(y) * d_3_0 + a_11(y) * f_1_0 + a_12(y) * f_3_0;
    let der_psi_2_0 = a_9(y) * d_2_0 + a_10(y) * d_4_0 + a_11(y) * f_2_0 + a_12(y) * f_4_0;
    let der_psi_3_0 = a_13(y) * d_1_0 + a_14(y) * d_3_0 + a_15(y) * f_1_0 + a_16(y) * f_3_0;
    let der_psi_4_0 = a_13(y) * d_2_0 + a_14(y) * d_4_0 + a_15(y) * f_2_0 + a_16(y) * f_4_0;

    let der_psi_1_1 = a_9(y) * d_1_1 + a_10(y) * d_3_1 + a_11(y) * f_1_1 + a_12(y) * f_3_1;
    let der_psi_2_1 = a_9(y) * d_2_1 + a_10(y) * d_4_1 + a_11(y) * f_2_1 + a_12(y) * f_4_1;
    let der_psi_3_1 = a_13(y) * d_1_1 + a_14(y) * d_3_1 + a_15(y) * f_1_1 + a_16(y) * f_3_1;
    let der_psi_4_1 = a_13(y) * d_2_1 + a_14(y) * d_4_1 + a_15(y) * f_2_1 + a_16(y) * f_4_1;

    (
        (der_psi_1_0, der_psi_2_0, der_psi_3_0, der_psi_4_0),
        (der_psi_1_1, der_psi_2_1, der_psi_3_1, der_psi_4_1),
    )
}

fn greens(
    y: f64,
    xi: f64,
    b: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
) -> (f64, f64, f64, f64) {
    let (
        (y_psi_1_0, y_psi_2_0, y_psi_3_0, y_psi_4_0),
        (y_psi_1_1, y_psi_2_1, y_psi_3_1, y_psi_4_1),
    ) = psi(y, b, alpha, mu_0, g, lambda);
    let (
        (xi_psi_1_0, xi_psi_2_0, xi_psi_3_0, xi_psi_4_0),
        (xi_psi_1_1, xi_psi_2_1, xi_psi_3_1, xi_psi_4_1),
    ) = psi(xi, b, alpha, mu_0, g, lambda);

    let g1 = if y < xi {
        y_psi_1_0 * xi_psi_1_1 + y_psi_2_0 * xi_psi_3_1
    } else {
        y_psi_1_1 * xi_psi_1_0 + y_psi_2_1 * xi_psi_3_0
    };
    let g2 = if y < xi {
        y_psi_1_0 * xi_psi_2_1 + y_psi_2_0 * xi_psi_4_1
    } else {
        y_psi_1_1 * xi_psi_2_0 + y_psi_2_1 * xi_psi_4_0
    };
    let g3 = if y < xi {
        y_psi_3_0 * xi_psi_1_1 + y_psi_4_0 * xi_psi_3_1
    } else {
        y_psi_3_1 * xi_psi_1_0 + y_psi_4_1 * xi_psi_3_0
    };
    let g4 = if y < xi {
        y_psi_3_0 * xi_psi_2_1 + y_psi_4_0 * xi_psi_4_1
    } else {
        y_psi_3_1 * xi_psi_2_0 + y_psi_4_1 * xi_psi_4_0
    };
    (g1, g2, g3, g4)
}

fn greens_dy(
    y: f64,
    xi: f64,
    b: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
) -> (f64, f64, f64, f64) {
    let (
        (y_psi_1_0, y_psi_2_0, y_psi_3_0, y_psi_4_0),
        (y_psi_1_1, y_psi_2_1, y_psi_3_1, y_psi_4_1),
    ) = der_psi(y, b, alpha, mu_0, g, lambda);
    let (
        (xi_psi_1_0, xi_psi_2_0, xi_psi_3_0, xi_psi_4_0),
        (xi_psi_1_1, xi_psi_2_1, xi_psi_3_1, xi_psi_4_1),
    ) = psi(xi, b, alpha, mu_0, g, lambda);

    let g1 = if y < xi {
        y_psi_1_0 * xi_psi_1_1 + y_psi_2_0 * xi_psi_3_1
    } else {
        y_psi_1_1 * xi_psi_1_0 + y_psi_2_1 * xi_psi_3_0
    };
    let g2 = if y < xi {
        y_psi_1_0 * xi_psi_2_1 + y_psi_2_0 * xi_psi_4_1
    } else {
        y_psi_1_1 * xi_psi_2_0 + y_psi_2_1 * xi_psi_4_0
    };
    let g3 = if y < xi {
        y_psi_3_0 * xi_psi_1_1 + y_psi_4_0 * xi_psi_3_1
    } else {
        y_psi_3_1 * xi_psi_1_0 + y_psi_4_1 * xi_psi_3_0
    };
    let g4 = if y < xi {
        y_psi_3_0 * xi_psi_2_1 + y_psi_4_0 * xi_psi_4_1
    } else {
        y_psi_3_1 * xi_psi_2_0 + y_psi_4_1 * xi_psi_4_0
    };
    (g1, g2, g3, g4)
}

fn a_1<F: Fn(f64) -> f64 + Send + Sync>(
    y: f64,
    a: f64,
    b: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    const N: i32 = 10;

    let sum = (1..N)
        .into_par_iter()
        .map(|n| {
            let alpha = alpha_calc(a, n as f64);
            let ((_, psi_2_0, _, _), _) = psi(y, b, alpha, mu_0, g, lambda);
            let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
            let pn =
                definite_integral(0.0, a, 50, eps, &|x| load_function(x) * f64::cos(alpha * x));

            sign * psi_2_0 * pn
        })
        .sum::<f64>();

    2.0 * sum / (1.0 + mu_0) / a
}

fn a_2(y: f64, xi: f64, a: f64, b: f64, mu_0: f64, g: f64, lambda: f64) -> f64 {
    const N: i32 = 20;

    let sum1 = (1..N)
        .into_par_iter()
        .map(|n| {
            let alpha = alpha_calc(a, n as f64);
            let (g1, _, _, _) = greens(y, xi, b, alpha, mu_0, g, lambda);
            g1
        })
        .sum::<f64>();

    let coef2 = 2.0 * a / PI;
    let sum2 = (1..N)
        .into_par_iter()
        .map(|n| {
            let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
            let a1 = 1.0 / (2.0 * n as f64 + 1.0);
            let a2 = f64::sin((2.0 * n as f64 + 1.0) * PI / 2.0);
            let a3 = f64::exp(-(2.0 * n as f64 + 1.0) * PI / 2.0 / a * (2.0 * b - y - xi));
            sign * a1 * a2 * a3
        })
        .sum::<f64>();

    let sum3 = a / 2.0 / PI * f64::ln(f64::cosh(PI / 2.0 / a * (2.0 * b - y - xi)) + 1.0);

    2.0 * (sum1 + coef2 * sum2 + sum3) / (1.0 + mu_0) / a
}

fn f_m<F: Fn(f64) -> f64 + Send + Sync>(
    m: usize,
    a: f64,
    b: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let f = |y| {
        let a_1 = a_1(y, a, b, mu_0, g, lambda, load_function, eps);
        let cheb = chebyshev(y, m);
        cheb * a_1
    };
    sqrt_gauss_integral(10, eps, &f)
}

fn g_k_m(k: usize, m: usize, a: f64, b: f64, mu_0: f64, g: f64, lambda: f64, eps: f64) -> f64 {
    let f = |y| {
        let cheb = chebyshev(y, m);
        let f = |xi| {
            let cheb = chebyshev(xi, k);
            let a_2 = a_2(y, xi, a, b, mu_0, g, lambda);
            a_2 * cheb
        };
        cheb * sqrt_gauss_integral(10, eps, &f)
    };
    let res = sqrt_gauss_integral(10, eps, &f) / PI;
    res
}

const N: usize = 5;
fn phi<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> Matrix<N, 1> {
    let left = (0..N)
        .into_par_iter()
        .map(|m| {
            let fm = f_m(m, a, b, mu_0, g, lambda, &load_function, eps);
            fm
        })
        .collect();
    let left = Matrix::<N, 1>::from_vec(left);

    let right = (0..N)
        .into_par_iter()
        .map(|m| {
            (0..N)
                .into_par_iter()
                .map(|k| {
                    let gkm = g_k_m(k, m, a, b, mu_0, g, lambda, eps);
                    if m == k {
                        let v_m = if m >= 1 { 1.0 / m as f64 } else { f64::ln(2.0) };
                        gkm + v_m * PI / 2.0
                    } else {
                        gkm
                    }
                })
                .collect::<Vec<_>>()
        })
        .flatten()
        .collect();
    let right = Matrix::<N, N>::from_vec(right);

    right.qr().solve(&left).unwrap()
}

pub fn unknown_function<F: Fn(f64) -> f64 + Send + Sync>(
    xi: f64,
    a: f64,
    b: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let phi = phi(a, b, mu_0, g, lambda, &load_function, eps);

    // let sqrt = f64::sqrt(1.0 - xi * xi);
    let sum: f64 = (0..N)
        .into_par_iter()
        .map(|k| phi.get(k).unwrap_or(&0.0) * chebyshev(xi, k))
        .sum();

    sum
}

pub fn function_un<F: Fn(f64) -> f64 + Send + Sync>(
    y: f64,
    n: usize,
    a: f64,
    b: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let alpha = alpha_calc(a, n as f64);
    let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
    let pn = definite_integral_limit(0.0, a, 10, &|x| load_function(x) * f64::cos(alpha * x));
    let ((_, psi_2_0, _, _), _) = psi(y, b, alpha, mu_0, g, lambda);

    let f = |xi: f64| {
        let xi = b / 2.0 * xi + b / 2.0;

        let unknown_fn = unknown_function(xi, a, b, mu_0, g, lambda, &load_function, eps);
        let (g1, _, _, _) = greens(y, xi, b, alpha, mu_0, g, lambda);
        g1 * unknown_fn
    };
    let int_val = sqrt_gauss_integral_finit(15, &f);

    let coef = sign * (1.0 + mu_0) * b / 2.0;

    coef * int_val - psi_2_0 * pn
}

pub fn function_vn<F: Fn(f64) -> f64 + Send + Sync>(
    y: f64,
    n: usize,
    a: f64,
    b: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let alpha = alpha_calc(a, n as f64);
    let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
    let pn = definite_integral_limit(0.0, a, 10, &|x| load_function(x) * f64::cos(alpha * x));
    let ((_, _, _, psi_4_0), _) = psi(y, b, alpha, mu_0, g, lambda);

    let f = |xi: f64| {
        let xi = b / 2.0 * xi + b / 2.0;

        let unknown_fn = unknown_function(xi, a, b, mu_0, g, lambda, &load_function, eps);
        let (_, _, g3, _) = greens(y, xi, b, alpha, mu_0, g, lambda);
        g3 * unknown_fn
    };
    let int_val = sqrt_gauss_integral_finit(15, &f);

    let coef = sign * (1.0 + mu_0) * b / 2.0;

    coef * int_val - psi_4_0 * pn
}

fn function_der_vn<F: Fn(f64) -> f64 + Send + Sync>(
    y: f64,
    n: usize,
    a: f64,
    b: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let alpha = alpha_calc(a, n as f64);
    let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
    let pn = definite_integral_limit(0.0, a, 10, &|x| load_function(x) * f64::cos(alpha * x));
    let ((_, _, _, psi_4_0), _) = psi(y, b, alpha, mu_0, g, lambda);

    let f = |xi: f64| {
        let xi = b / 2.0 * xi + b / 2.0;

        let unknown_fn = unknown_function(xi, a, b, mu_0, g, lambda, &load_function, eps);
        let (_, _, g3, _) = greens_dy(y, xi, b, alpha, mu_0, g, lambda);
        g3 * unknown_fn
    };
    let int_val = sqrt_gauss_integral_finit(15, &f);

    let coef = sign * (1.0 + mu_0) * b / 2.0;

    coef * int_val - psi_4_0 * pn
}

pub fn function_u<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    x: f64,
    y: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let n = 3;
    let start = 1;
    let f = |n| {
        let alpha = alpha_calc(a, n as f64);
        2.0 * function_un(y, n, a, b, mu_0, g, lambda, load_function, eps) * f64::sin(alpha * x) / a
    };

    sum_calc_finit(&f, start, n)
}

pub fn function_u_dx<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    x: f64,
    y: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let n = 3;
    let start = 1;
    let f = |n| {
        let alpha = alpha_calc(a, n as f64);
        2.0 * alpha
            * function_un(y, n, a, b, mu_0, g, lambda, load_function, eps)
            * f64::cos(alpha * x)
            / a
    };

    sum_calc_finit(&f, start, n)
}

pub fn function_v<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    x: f64,
    y: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let n = 3;
    let start = 1;
    let f = |n| {
        let alpha = alpha_calc(a, n as f64);
        2.0 * function_vn(y, n, a, b, mu_0, g, lambda, load_function, eps) * f64::cos(alpha * x) / a
    };

    sum_calc_finit(&f, start, n)
}

pub fn function_v_dy<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    x: f64,
    y: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let n = 3;
    let start = 1;
    let f = |n| {
        let alpha = alpha_calc(a, n as f64);
        2.0 * function_der_vn(y, n, a, b, mu_0, g, lambda, load_function, eps) * f64::cos(alpha * x)
            / a
    };

    sum_calc_finit(&f, start, n)
}

fn function_sigma_x<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    x: f64,
    y: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let d_ux = function_u_dx(a, b, x, y, mu_0, g, lambda, load_function, eps);
    let d_vy = function_v_dy(a, b, x, y, mu_0, g, lambda, load_function, eps);

    2_f64 * g * d_ux + lambda * d_vy + lambda * d_ux
}

fn function_sigma_y<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    x: f64,
    y: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let d_ux = function_u_dx(a, b, x, y, mu_0, g, lambda, load_function, eps);
    let d_vy = function_v_dy(a, b, x, y, mu_0, g, lambda, load_function, eps);

    (2.0 * g + lambda) * d_vy + lambda * d_ux
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::{g, lambda, mu_0};
    use nalgebra::{Matrix2, RowVector2};

    #[test]
    fn coefficients_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;
        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);
        let alpha = PI / a * 2.0;

        let (b_1, b_2, b_3, b_4, b_5, b_6, b_7, b_8, b_9, b_10, b_11, b_12, b_13, b_14, b_15, b_16) =
            b_values(b, alpha, mu_0, g, lambda);

        let eq_1 =
            |x_1: f64, x_2: f64, x_3: f64, x_4: f64| b_1 * x_1 + b_2 * x_2 + b_3 * x_3 + b_4 * x_4;
        let eq_2 =
            |x_1: f64, x_2: f64, x_3: f64, x_4: f64| b_5 * x_1 + b_6 * x_2 + b_7 * x_3 + b_8 * x_4;
        let eq_3 = |x_1: f64, x_2: f64, x_3: f64, x_4: f64| {
            b_9 * x_1 + b_10 * x_2 + b_11 * x_3 + b_12 * x_4
        };
        let eq_4 = |x_1: f64, x_2: f64, x_3: f64, x_4: f64| {
            b_13 * x_1 + b_14 * x_2 + b_15 * x_3 + b_16 * x_4
        };

        let (
            ((d_1_0, d_2_0, d_3_0, d_4_0), (f_1_0, f_2_0, f_3_0, f_4_0)),
            ((d_1_1, d_2_1, d_3_1, d_4_1), (f_1_1, f_2_1, f_3_1, f_4_1)),
        ) = coefficients(b, alpha, mu_0, g, lambda);

        println!("system 1");
        println!("{}", eq_1(d_1_0, d_3_0, f_1_0, f_3_0));
        println!("{}", eq_2(d_1_0, d_3_0, f_1_0, f_3_0));
        println!("{}", eq_3(d_1_0, d_3_0, f_1_0, f_3_0));
        println!("{}", eq_4(d_1_0, d_3_0, f_1_0, f_3_0));

        println!("system 2");
        println!("{}", eq_1(d_2_0, d_4_0, f_2_0, f_4_0));
        println!("{}", eq_2(d_2_0, d_4_0, f_2_0, f_4_0));
        println!("{}", eq_3(d_2_0, d_4_0, f_2_0, f_4_0));
        println!("{}", eq_4(d_2_0, d_4_0, f_2_0, f_4_0));

        println!("system 3");
        println!("{}", eq_1(d_1_1, d_3_1, f_1_1, f_3_1));
        println!("{}", eq_2(d_1_1, d_3_1, f_1_1, f_3_1));
        println!("{}", eq_3(d_1_1, d_3_1, f_1_1, f_3_1));
        println!("{}", eq_4(d_1_1, d_3_1, f_1_1, f_3_1));

        println!("system 4");
        println!("{}", eq_1(d_2_1, d_4_1, f_2_1, f_4_1));
        println!("{}", eq_2(d_2_1, d_4_1, f_2_1, f_4_1));
        println!("{}", eq_3(d_2_1, d_4_1, f_2_1, f_4_1));
        println!("{}", eq_4(d_2_1, d_4_1, f_2_1, f_4_1));
    }

    #[test]
    fn psi_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;
        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);
        let alpha = PI / a * 2.0;

        let (
            (b_psi_1_0, b_psi_2_0, b_psi_3_0, b_psi_4_0),
            (b_psi_1_1, b_psi_2_1, b_psi_3_1, b_psi_4_1),
        ) = psi(b, b, alpha, mu_0, g, lambda);
        let (
            (b_der_psi_1_0, b_der_psi_2_0, b_der_psi_3_0, b_der_psi_4_0),
            (b_der_psi_1_1, b_der_psi_2_1, b_der_psi_3_1, b_der_psi_4_1),
        ) = der_psi(b, b, alpha, mu_0, g, lambda);
        let (
            (zero_psi_1_0, zero_psi_2_0, zero_psi_3_0, zero_psi_4_0),
            (zero_psi_1_1, zero_psi_2_1, zero_psi_3_1, zero_psi_4_1),
        ) = psi(0.0, b, alpha, mu_0, g, lambda);
        let (
            (zero_der_psi_1_0, zero_der_psi_2_0, zero_der_psi_3_0, zero_der_psi_4_0),
            (zero_der_psi_1_1, zero_der_psi_2_1, zero_der_psi_3_1, zero_der_psi_4_1),
        ) = der_psi(0.0, b, alpha, mu_0, g, lambda);

        let a_1 = Matrix2::from_rows(&[
            RowVector2::new(1.0, 0.0),
            RowVector2::new(0.0, 2.0 * g + lambda),
        ]);
        let a_2 = Matrix2::from_rows(&[
            RowVector2::new(0.0, -alpha),
            RowVector2::new(alpha * lambda, 0.0),
        ]);

        let b_der_psi_0 = Matrix2::from_rows(&[
            RowVector2::new(b_der_psi_1_0, b_der_psi_2_0),
            RowVector2::new(b_der_psi_3_0, b_der_psi_4_0),
        ]);
        let b_psi_0 = Matrix2::from_rows(&[
            RowVector2::new(b_psi_1_0, b_psi_2_0),
            RowVector2::new(b_psi_3_0, b_psi_4_0),
        ]);
        let res = a_1 * b_der_psi_0 + a_2 * b_psi_0;
        println!("U_0[Psi_0]: {}", res);

        let b_der_psi_0 = Matrix2::from_rows(&[
            RowVector2::new(b_der_psi_1_1, b_der_psi_2_1),
            RowVector2::new(b_der_psi_3_1, b_der_psi_4_1),
        ]);
        let b_psi_0 = Matrix2::from_rows(&[
            RowVector2::new(b_psi_1_1, b_psi_2_1),
            RowVector2::new(b_psi_3_1, b_psi_4_1),
        ]);
        let res = a_1 * b_der_psi_0 + a_2 * b_psi_0;
        println!("U_0[Psi_1]: {}", res);

        let a_1 = Matrix2::from_rows(&[RowVector2::new(1.0, 0.0), RowVector2::new(0.0, 0.0)]);
        let a_2 = Matrix2::from_rows(&[RowVector2::new(0.0, -alpha), RowVector2::new(0.0, 1.0)]);

        let zero_der_psi_0 = Matrix2::from_rows(&[
            RowVector2::new(zero_der_psi_1_0, zero_der_psi_2_0),
            RowVector2::new(zero_der_psi_3_0, zero_der_psi_4_0),
        ]);
        let zero_psi_0 = Matrix2::from_rows(&[
            RowVector2::new(zero_psi_1_0, zero_psi_2_0),
            RowVector2::new(zero_psi_3_0, zero_psi_4_0),
        ]);
        let res = a_1 * zero_der_psi_0 + a_2 * zero_psi_0;
        println!("U_1[Psi_0]: {}", res);

        let zero_der_psi_0 = Matrix2::from_rows(&[
            RowVector2::new(zero_der_psi_1_1, zero_der_psi_2_1),
            RowVector2::new(zero_der_psi_3_1, zero_der_psi_4_1),
        ]);
        let zero_psi_0 = Matrix2::from_rows(&[
            RowVector2::new(zero_psi_1_1, zero_psi_2_1),
            RowVector2::new(zero_psi_3_1, zero_psi_4_1),
        ]);
        let res = a_1 * zero_der_psi_0 + a_2 * zero_psi_0;
        println!("U_1[Psi_1]: {}", res);
    }

    #[test]
    fn psi_lim_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;
        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);

        for n in 1..50 {
            let alpha = PI / a * (n as f64 - 0.5);

            let ((psi_1_0, psi_2_0, psi_3_0, psi_4_0), (psi_1_1, psi_2_1, psi_3_1, psi_4_1)) =
                psi(b, b, alpha, mu_0, g, lambda);

            println!("alpha: {alpha}");
            println!(
                "psi_1_0: {psi_1_0}, psi_2_0: {psi_2_0}, psi_3_0: {psi_3_0}, psi_4_0: {psi_4_0}"
            );
            println!(
                "psi_1_1: {psi_1_1}, psi_2_1: {psi_2_1}, psi_3_1: {psi_3_1}, psi_4_1: {psi_4_1}"
            );
        }
    }

    #[test]
    fn a1_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;
        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);
        let load_function = |x| x * x;
        let eps = 0.01;

        let y = 5.5;

        let a1 = a_1(y, a, b, mu_0, g, lambda, &load_function, eps);
        println!("a1: {}", a1);
    }

    #[test]
    fn a2_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;
        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);

        let y = 5.5;
        let xi = 6.5;

        let a2 = a_2(y, xi, a, b, mu_0, g, lambda);
        println!("a2: {}", a2);
    }

    #[test]
    fn fm_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;
        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);
        let load_function = |x| x * x;
        let eps = 0.1;

        for m in 0..1000 {
            let f_m = f_m(m, a, b, mu_0, g, lambda, &load_function, eps);
            println!("m:{m}, {f_m}");
        }
    }

    #[test]
    fn gkm_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;
        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);
        let eps = 0.1;

        for k in 0..10 {
            for m in 0..10 {
                let g_k_m = g_k_m(k, m, a, b, mu_0, g, lambda, eps);
                println!("{g_k_m}");
            }
        }
    }

    #[test]
    fn phi_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;
        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);
        let load_function = |x| x * x;
        let eps = 0.1;

        let phi = phi(a, b, mu_0, g, lambda, &load_function, eps);
        println!("phi: {}", phi);
    }

    #[test]
    fn unknown_function_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;
        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);
        let load_function = |x| x * x;
        let eps = 0.1;

        let xi = 1.0;

        let unknown_function = unknown_function(xi, a, b, mu_0, g, lambda, &load_function, eps);
        println!("unknown_function: {}", unknown_function);
    }

    #[test]
    fn random_test() {
        let a = 10.0;
        let b = 15.0;
        // steel
        let puasson_coef = 0.25;
        let young_modulus = 200.0;

        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);

        let load_function = |x| x * x;
        let eps = 0.1;

        let y = b;
        let x = 0.5;

        let res = function_un(y, 1, a, b, mu_0, g, lambda, &load_function, eps);

        println!("res: {}", res);
    }
}
