use std::f64::consts::PI;

use crate::integration::definite_integral;

fn a_coefficients(omega: f64, c1: f64, c2: f64, alpha: f64, mu_0: f64) -> (f64, f64) {
    let x1 = 2_f64 * alpha * alpha * c1 * c1 * c2 * c2 * mu_0
        + 2_f64 * alpha * alpha * c1 * c1 * c2 * c2
        - c1 * c1 * omega * omega
        - c2 * c2 * mu_0 * omega * omega
        - c2 * c2 * omega * omega;

    let x2 = 4_f64 * alpha * alpha * c1 * c1 * c1 * c1 * c2 * c2 * mu_0 * mu_0
        + 4_f64 * alpha * alpha * c1 * c1 * c1 * c1 * c2 * c2 * mu_0
        - 4_f64 * alpha * alpha * c1 * c1 * c2 * c2 * c2 * c2 * mu_0 * mu_0
        - 4_f64 * alpha * alpha * c1 * c1 * c2 * c2 * c2 * c2 * mu_0
        + c1 * c1 * c1 * c1 * omega * omega
        - 2_f64 * c1 * c1 * c2 * c2 * mu_0 * omega * omega
        - 2_f64 * c1 * c1 * c2 * c2 * omega * omega
        + c2 * c2 * c2 * c2 * mu_0 * mu_0 * omega * omega
        + 2_f64 * c2 * c2 * c2 * c2 * omega * omega * mu_0
        + c2 * c2 * c2 * c2 * omega * omega;

    let x3 = 2_f64 * c1 * c1 * c2 * c2 * mu_0 + 2_f64 * c1 * c1 * c2 * c2;

    let a1 = -f64::sqrt((x1 - omega * f64::sqrt(x2)) / x3);
    let a3 = -f64::sqrt((x1 + omega * f64::sqrt(x2)) / x3);

    (a1, a3)
}

fn c_coefficients<F: Fn(f64) -> f64>(
    a: f64,
    b: f64,
    a1: f64,
    a3: f64,
    omega: f64,
    c1: f64,
    c2: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> (f64, f64, f64, f64) {
    let pn = definite_integral(0_f64, a, 100, eps, &|x| {
        load_function(x) * f64::cos(alpha * x)
    });
    let coef = pn * (1_f64 + mu_0);

    let e1 = f64::exp(a1 * b) + f64::exp(-a1 * b);
    let e2 = f64::exp(a1 * b) - f64::exp(-a1 * b);
    let e3 = f64::exp(a3 * b) + f64::exp(-a3 * b);
    let e4 = f64::exp(a3 * b) - f64::exp(-a3 * b);

    let z1 = 2_f64 * a1 * (a1 * a1 - a3 * a3);
    let z2 = 2_f64 * a3 * (a3 * a3 - a1 * a1);

    let x1 = a1 * alpha * mu_0;
    let x2 = a3 * alpha * mu_0;
    let _x3 = a1 * a1 + a1 * a1 * mu_0 - alpha * alpha + omega * omega / c2 / c2;
    let x4 = a3 * a3 + a3 * a3 * mu_0 - alpha * alpha + omega * omega / c2 / c2;
    let x5 = a1 * a1 - alpha * alpha - alpha * alpha * mu_0 + omega * omega / c1 / c1;
    let x6 = a3 * a3 - alpha * alpha - alpha * alpha * mu_0 + omega * omega / c1 / c1;

    let d2 = e2 * (a1 * x1 - alpha * x5) / z1;
    let d4 = e4 * (a3 * x4 - alpha * x6) / z2;
    let d6 = e1 * (a1 * x5 * (2_f64 * g + lambda) + alpha * lambda * x1) / z1;
    let d8 = e3 * (a3 * x6 * (2_f64 * g + lambda) + alpha * lambda * x2) / z2;

    let c1 = 0_f64;
    let c2 = coef * d4 / (d8 * d2 - d4 * d6);
    let c3 = 0_f64;
    let c4 = -coef * d2 / (d8 * d2 - d4 * d6);

    (c1, c2, c3, c4)
}

fn function_un<F: Fn(f64) -> f64>(
    y: f64,
    a: f64,
    b: f64,
    omega: f64,
    c1: f64,
    c2: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let (a1, a3) = a_coefficients(omega, c1, c2, alpha, mu_0);

    let e1 = f64::exp(a1 * y);
    let e2 = f64::exp(-a1 * y);
    let e3 = f64::exp(a3 * y);
    let e4 = f64::exp(-a3 * y);

    let (c1, c2, c3, c4) = c_coefficients(
        a,
        b,
        a1,
        a3,
        omega,
        c1,
        c2,
        alpha,
        mu_0,
        g,
        lambda,
        load_function,
        eps,
    );

    let res = c1 * (a1 * a1 + a1 * a1 * mu_0 - alpha * alpha + omega * omega / c2 / c2) * (e1 - e2)
        / (2_f64 * a1 * (a1 * a1 - a3 * a3))
        + c2 * a1 * alpha * mu_0 * (e1 + e2) / (2_f64 * a1 * (a1 * a1 - a3 * a3))
        + c3 * (a3 * a3 + a3 * a3 * mu_0 - alpha * alpha + omega * omega / c2 / c2) * (e3 - e4)
            / (2_f64 * a3 * (a3 * a3 - a1 * a1))
        + c4 * a3 * alpha * mu_0 * (e3 + e4) / (2_f64 * a3 * (a3 * a3 - a1 * a1));

    res
}

fn function_vn<F: Fn(f64) -> f64>(
    y: f64,
    a: f64,
    b: f64,
    omega: f64,
    c1: f64,
    c2: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let (a1, a3) = a_coefficients(omega, c1, c2, alpha, mu_0);

    let e1 = f64::exp(a1 * y);
    let e2 = f64::exp(-a1 * y);
    let e3 = f64::exp(a3 * y);
    let e4 = f64::exp(-a3 * y);

    let (c1, c2, c3, c4) = c_coefficients(
        a,
        b,
        a1,
        a3,
        omega,
        c1,
        c2,
        alpha,
        mu_0,
        g,
        lambda,
        load_function,
        eps,
    );

    let res = c1 * (-a1 * alpha * mu_0) * (e1 + e2) / (2_f64 * a1 * (a1 * a1 - a3 * a3))
        + c2 * (a1 * a1 - alpha * alpha - alpha * alpha * mu_0 + omega * omega / c1 / c1)
            * (e1 - e2)
            / (2_f64 * a1 * (a1 * a1 - a3 * a3))
        + c3 * (-a3 * alpha * mu_0) * (e3 + e4) / (2_f64 * a3 * (a3 * a3 - a1 * a1))
        + c4 * (a3 * a3 - alpha * alpha - alpha * alpha * mu_0 + omega * omega / c1 / c1)
            * (e3 - e4)
            / (2_f64 * a3 * (a3 * a3 - a1 * a1));

    res
}

pub fn function_u<F: Fn(f64) -> f64>(
    a: f64,
    b: f64,
    x: f64,
    y: f64,
    omega: f64,
    c1: f64,
    c2: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let mut n = 10;
    let mut result = 0_f64;
    let mut prev_result;

    if x != a && x != 0_f64 {
        for i in 1..n {
            let alpha = PI * i as f64 / a;
            result += 2_f64
                * function_un(
                    y,
                    a,
                    b,
                    omega,
                    c1,
                    c2,
                    alpha,
                    mu_0,
                    g,
                    lambda,
                    load_function,
                    eps,
                )
                * f64::sin(alpha * x)
                / a
                / (1_f64 + mu_0);
        }
    }

    loop {
        n *= 2;
        prev_result = result;
        result = 0_f64;

        if x != a && x != 0_f64 {
            for i in 1..n {
                let alpha = PI * i as f64 / a;
                result += 2_f64
                    * function_un(
                        y,
                        a,
                        b,
                        omega,
                        c1,
                        c2,
                        alpha,
                        mu_0,
                        g,
                        lambda,
                        load_function,
                        eps,
                    )
                    * f64::sin(alpha * x)
                    / a
                    / (1_f64 + mu_0);
            }
        }

        if f64::abs(result - prev_result) < eps {
            break;
        }
    }

    result
}
