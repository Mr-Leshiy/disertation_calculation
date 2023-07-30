use crate::{integration::definite_integral, utils::sum_calc};
use num::{complex::Complex64, Zero};
use std::f64::consts::PI;

fn det_roots(omega: f64, c1: f64, c2: f64, alpha: f64, mu_0: f64) -> (Complex64, Complex64) {
    let a1 = -2.0 * alpha * alpha + omega * omega / c2 / c2 - 2.0 * alpha * alpha * mu_0
        + omega * omega / c1 / c1
        + mu_0 * omega * omega / c1 / c1;
    let a2 = alpha * alpha * alpha * alpha - alpha * alpha * omega * omega / c2 / c2
        + alpha * alpha * alpha * alpha * mu_0
        - alpha * alpha * mu_0 * omega * omega / c2 / c2
        - alpha * alpha * omega * omega / c2 / c2
        + omega * omega * omega * omega / c1 / c1 / c2 / c2;

    let sqrt = f64::sqrt(a1 * a1 - 4.0 * (1.0 + mu_0) * a2);
    let s1 = -a1 + sqrt;
    let s2 = -a1 - sqrt;
    let coef = 2.0 * (1.0 + mu_0);
    let s1 = if s1 >= 0.0 {
        Complex64::new(f64::sqrt(s1 / coef), 0.0)
    } else {
        Complex64::new(0.0, f64::sqrt(-s1 / coef))
    };
    let s2 = if s2 >= 0.0 {
        Complex64::new(f64::sqrt(s2 / coef), 0.0)
    } else {
        Complex64::new(0.0, f64::sqrt(-s2 / coef))
    };

    (s1, s2)
}

fn c_coefficients<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    s1: Complex64,
    s2: Complex64,
    omega: f64,
    c1: f64,
    c2: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> (Complex64, Complex64, Complex64, Complex64) {
    let pn = definite_integral(0.0, a, 100, eps, &|x| {
        load_function(x) * f64::cos(alpha * x)
    });
    let coef = pn * (1.0 + mu_0);

    let e1 = Complex64::exp(s1 * b) + Complex64::exp(-s1 * b);
    let e2 = Complex64::exp(s1 * b) - Complex64::exp(-s1 * b);
    let e3 = Complex64::exp(s2 * b) + Complex64::exp(-s2 * b);
    let e4 = Complex64::exp(s2 * b) - Complex64::exp(-s2 * b);

    let z1 = 2.0 * s1 * (s1 * s1 - s2 * s2);
    let z2 = 2.0 * s2 * (s2 * s2 - s1 * s1);

    let x1 = s1 * alpha * mu_0;
    let x2 = s2 * alpha * mu_0;
    let _x3 = s1 * s1 + s1 * s1 * mu_0 - alpha * alpha + omega * omega / c2 / c2;
    let x4 = s2 * s2 + s2 * s2 * mu_0 - alpha * alpha + omega * omega / c2 / c2;
    let x5 = s1 * s1 - alpha * alpha - alpha * alpha * mu_0 + omega * omega / c1 / c1;
    let x6 = s2 * s2 - alpha * alpha - alpha * alpha * mu_0 + omega * omega / c1 / c1;

    let d2 = e2 * (s1 * x1 - alpha * x5) / z1;
    let d4 = e4 * (s2 * x4 - alpha * x6) / z2;
    let d6 = e1 * (s1 * x5 * (2.0 * g + lambda) + alpha * lambda * x1) / z1;
    let d8 = e3 * (s2 * x6 * (2.0 * g + lambda) + alpha * lambda * x2) / z2;

    let c1 = Complex64::zero();
    let c2 = coef * d4 / (d8 * d2 - d4 * d6);
    let c3 = Complex64::zero();
    let c4 = -coef * d2 / (d8 * d2 - d4 * d6);

    (c1, c2, c3, c4)
}

fn function_un<F: Fn(f64) -> f64 + Send + Sync>(
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
) -> Complex64 {
    let (a1, a3) = det_roots(omega, c1, c2, alpha, mu_0);

    let e1 = Complex64::exp(a1 * y);
    let e2 = Complex64::exp(-a1 * y);
    let e3 = Complex64::exp(a3 * y);
    let e4 = Complex64::exp(-a3 * y);

    let (c_1, c_2, c_3, c_4) = c_coefficients(
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

    let res = c_1
        * (a1 * a1 + a1 * a1 * mu_0 - alpha * alpha + omega * omega / c2 / c2)
        * (e1 - e2)
        / (2.0 * a1 * (a1 * a1 - a3 * a3))
        + c_2 * a1 * alpha * mu_0 * (e1 + e2) / (2.0 * a1 * (a1 * a1 - a3 * a3))
        + c_3 * (a3 * a3 + a3 * a3 * mu_0 - alpha * alpha + omega * omega / c2 / c2) * (e3 - e4)
            / (2.0 * a3 * (a3 * a3 - a1 * a1))
        + c_4 * a3 * alpha * mu_0 * (e3 + e4) / (2.0 * a3 * (a3 * a3 - a1 * a1));

    res
}

fn function_vn<F: Fn(f64) -> f64 + Send + Sync>(
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
) -> Complex64 {
    let (s1, s2) = det_roots(omega, c1, c2, alpha, mu_0);

    let e1 = Complex64::exp(s1 * y);
    let e2 = Complex64::exp(-s1 * y);
    let e3 = Complex64::exp(s2 * y);
    let e4 = Complex64::exp(-s2 * y);

    let (c_1, c_2, c_3, c_4) = c_coefficients(
        a,
        b,
        s1,
        s2,
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

    let res = c_1 * (-s1 * alpha * mu_0) * (e1 + e2) / (2.0 * s1 * (s1 * s1 - s2 * s2))
        + c_2
            * (s1 * s1 - alpha * alpha - alpha * alpha * mu_0 + omega * omega / c1 / c1)
            * (e1 - e2)
            / (2.0 * s1 * (s1 * s1 - s2 * s2))
        + c_3 * (-s2 * alpha * mu_0) * (e3 + e4) / (2.0 * s2 * (s2 * s2 - s1 * s1))
        + c_4
            * (s2 * s2 - alpha * alpha - alpha * alpha * mu_0 + omega * omega / c1 / c1)
            * (e3 - e4)
            / (2.0 * s2 * (s2 * s2 - s1 * s1));

    res
}

fn function_derivative_vn<F: Fn(f64) -> f64 + Send + Sync>(
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
) -> Complex64 {
    let (s1, s2) = det_roots(omega, c1, c2, alpha, mu_0);

    let e1 = Complex64::exp(s1 * y);
    let e2 = Complex64::exp(-s1 * y);
    let e3 = Complex64::exp(s2 * y);
    let e4 = Complex64::exp(-s2 * y);

    let (c_1, c_2, c_3, c_4) = c_coefficients(
        a,
        b,
        s1,
        s2,
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

    let res = c_1 * (-s1 * alpha * mu_0) * s1 * (e1 - e2) / (2_f64 * s1 * (s1 * s1 - s2 * s2))
        + c_2
            * (s1 * s1 - alpha * alpha - alpha * alpha * mu_0 + omega * omega / c1 / c1)
            * s1
            * (e1 + e2)
            / (2_f64 * s1 * (s1 * s1 - s2 * s2))
        + c_3 * (-s2 * alpha * mu_0) * s2 * (e3 - e4) / (2_f64 * s2 * (s2 * s2 - s1 * s1))
        + c_4
            * (s2 * s2 - alpha * alpha - alpha * alpha * mu_0 + omega * omega / c1 / c1)
            * s2
            * (e3 + e4)
            / (2_f64 * s2 * (s2 * s2 - s1 * s1));

    res
}

fn function_u<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    x: f64,
    y: f64,
    t: f64,
    omega: f64,
    c1: f64,
    c2: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let n = 10;
    let start = 1;

    // actually we should have a e^(i*omega*t), but for simplicity we are taking only real part
    let cos_t = f64::cos(omega * t);

    let f = |i| {
        if x != a && x != 0_f64 {
            let alpha = PI * i as f64 / a;
            let res =
                2.0 * function_un(
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
                ) * f64::sin(alpha * x)
                    / a;
            res.re
        } else {
            0.0
        }
    };

    sum_calc(&f, eps, start, n) * cos_t
}

fn function_derivative_u_x<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    x: f64,
    y: f64,
    t: f64,
    omega: f64,
    c1: f64,
    c2: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let n = 10;
    let start = 1;

    // actually we should have a e^(i*omega*t), but for simplicity we are taking only real part
    let cos_t = f64::cos(omega * t);

    let f = |i| {
        if x != a && x != 0_f64 {
            let alpha = PI * i as f64 / a;
            let res = 2.0
                * alpha
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
                * f64::cos(alpha * x)
                / a;
            res.re
        } else {
            0.0
        }
    };

    sum_calc(&f, eps, start, n) * cos_t
}

fn function_v<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    x: f64,
    y: f64,
    t: f64,
    omega: f64,
    c1: f64,
    c2: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let p0 = definite_integral(0_f64, a, 100, eps, load_function);
    let v0 = -p0 * f64::sin(y * f64::sqrt(omega * omega / (c2 * c2 * (1_f64 + mu_0))))
        / ((2_f64 * g + lambda)
            * f64::sqrt(omega * omega / (c2 * c2 * (1_f64 + mu_0)))
            * f64::sin(b * f64::sqrt(omega * omega / (c2 * c2 * (1_f64 + mu_0)))));
    let v0 = v0 / a;
    let n = 10;
    let start = 1;

    // actually we should have a e^(i*omega*t), but for simplicity we are taking only real part
    let cos_t = f64::cos(omega * t);

    let f = |i| {
        let alpha = PI * i as f64 / a;
        let res =
            2.0 * function_vn(
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
            ) * f64::cos(alpha * x)
                / a;
        res.re
    };

    v0 * cos_t + sum_calc(&f, eps, start, n) * cos_t
}

fn function_derivative_v_y<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    x: f64,
    y: f64,
    t: f64,
    omega: f64,
    c1: f64,
    c2: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let p0 = definite_integral(0_f64, a, 100, eps, load_function);
    let v0 = -p0 * f64::sin(y * f64::sqrt(omega * omega / (c2 * c2 * (1_f64 + mu_0))))
        / ((2_f64 * g + lambda)
            * f64::sqrt(omega * omega / (c2 * c2 * (1_f64 + mu_0)))
            * f64::sin(b * f64::sqrt(omega * omega / (c2 * c2 * (1_f64 + mu_0)))));
    let v0 = v0 / a;
    let n = 10;
    let start = 1;

    // actually we should have a e^(i*omega*t), but for simplicity we are taking only real part
    let cos_t = f64::cos(omega * t);

    let f = |i| {
        let alpha = PI * i as f64 / a;
        let res =
            2.0 * function_derivative_vn(
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
            ) * f64::cos(alpha * x)
                / a;
        res.re
    };

    v0 * cos_t + sum_calc(&f, eps, start, n) * cos_t
}

fn function_sigma_x<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    x: f64,
    y: f64,
    t: f64,
    omega: f64,
    c1: f64,
    c2: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let d_ux = function_derivative_u_x(
        a,
        b,
        x,
        y,
        t,
        omega,
        c1,
        c2,
        mu_0,
        g,
        lambda,
        load_function,
        eps,
    );
    let d_vy = function_derivative_v_y(
        a,
        b,
        x,
        y,
        t,
        omega,
        c1,
        c2,
        mu_0,
        g,
        lambda,
        load_function,
        eps,
    );

    2_f64 * g * d_ux + lambda * d_vy + lambda * d_ux
}

fn function_sigma_y<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    x: f64,
    y: f64,
    t: f64,
    omega: f64,
    c1: f64,
    c2: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let d_ux = function_derivative_u_x(
        a,
        b,
        x,
        y,
        t,
        omega,
        c1,
        c2,
        mu_0,
        g,
        lambda,
        load_function,
        eps,
    );
    let d_vy = function_derivative_v_y(
        a,
        b,
        x,
        y,
        t,
        omega,
        c1,
        c2,
        mu_0,
        g,
        lambda,
        load_function,
        eps,
    );

    2_f64 * g * d_vy + lambda * d_vy + lambda * d_ux
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::{g, lambda, mu_0};

    #[test]
    fn det_roots_test() {
        let a: f64 = 10_f64;
        let omega = 0.75;
        let c1 = 10_f64;
        let c2 = 10_f64;
        let puasson_coef = 0.25;
        let eps = 0.001;

        let mu_0 = mu_0(puasson_coef);

        for i in 1..10 {
            let alpha = PI * i as f64 / a;

            let (s1, s2) = det_roots(omega, c1, c2, alpha, mu_0);
            let (s1, s2) = (s1.re, s2.re);

            let func = |s: f64| {
                (s * s - alpha * alpha - alpha * alpha * mu_0 + omega * omega / c1 / c1)
                    * (s * s + s * s * mu_0 - alpha * alpha + omega * omega / c2 / c2)
                    + s * s * alpha * alpha * mu_0 * mu_0
            };

            assert!(
                f64::abs(func(s1) - 0_f64) < eps,
                "result: {}, exp: {}",
                func(s1),
                0_f64
            );
            assert!(
                f64::abs(func(-s1) - 0_f64) < eps,
                "result: {}, exp: {}",
                func(-s1),
                0_f64
            );
            assert!(
                f64::abs(func(s2) - 0_f64) < eps,
                "result: {}, exp: {}",
                func(s2),
                0_f64
            );
            assert!(
                f64::abs(func(-s2) - 0_f64) < eps,
                "result: {}, exp: {}",
                func(-s2),
                0_f64
            );
        }
    }

    #[test]
    fn c_coefficients_test() {
        let load_function = |x| x * x;
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;
        let omega = 0.75;
        let c1 = 10_f64;
        let c2 = 10_f64;
        let eps = 0.001;

        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);

        for i in 1..10 {
            let alpha = PI * i as f64 / a;

            let (a1, a3) = det_roots(omega, c1, c2, alpha, mu_0);

            let pn = definite_integral(0_f64, a, 100, eps, &|x| {
                load_function(x) * f64::cos(alpha * x)
            });

            let e1 = Complex64::exp(a1 * b) + Complex64::exp(-a1 * b);
            let e2 = Complex64::exp(a1 * b) - Complex64::exp(-a1 * b);
            let e3 = Complex64::exp(a3 * b) + Complex64::exp(-a3 * b);
            let e4 = Complex64::exp(a3 * b) - Complex64::exp(-a3 * b);

            let z1 = 2_f64 * a1 * (a1 * a1 - a3 * a3);
            let z2 = 2_f64 * a3 * (a3 * a3 - a1 * a1);

            let x1 = a1 * alpha * mu_0;
            let x2 = a3 * alpha * mu_0;
            let x3 = a1 * a1 + a1 * a1 * mu_0 - alpha * alpha + omega * omega / c2 / c2;
            let x4 = a3 * a3 + a3 * a3 * mu_0 - alpha * alpha + omega * omega / c2 / c2;
            let x5 = a1 * a1 - alpha * alpha - alpha * alpha * mu_0 + omega * omega / c1 / c1;
            let x6 = a3 * a3 - alpha * alpha - alpha * alpha * mu_0 + omega * omega / c1 / c1;

            let d1 = e1 * (a1 * x3 + alpha * x1) / z1;
            let d2 = e2 * (a1 * x1 - alpha * x5) / z1;
            let d3 = e3 * (a3 * x4 + alpha * x2) / z2;
            let d4 = e4 * (a3 * x4 - alpha * x6) / z2;
            let d5 = e2 * (a1 * x3 - a1 * x1 * (2_f64 * g + lambda)) / z1;
            let d6 = e1 * (a1 * x5 * (2_f64 * g + lambda) + alpha * lambda * x1) / z1;
            let d7 = e4 * (alpha * lambda * x4 - a3 * x2 * (2_f64 * g + lambda)) / z2;
            let d8 = e3 * (a3 * x6 * (2_f64 * g + lambda) + alpha * lambda * x2) / z2;

            let (c_1, c_2, c_3, c_4) = c_coefficients(
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
                &load_function,
                eps,
            );

            let eq1 = c_1 * d1 + c_2 * d2 + c_3 * d3 + c_4 * d4;
            let eq2 = c_1 * d5 + c_2 * d6 + c_3 * d7 + c_4 * d8;
            let eq3 = c_1 * ((a1 * x3 + alpha * x1) / z1) + c_3 * ((alpha * x2 + a3 * x4) / z2);
            let eq4 = c_1 * (-x1 / z1) + c_3 * (-x2 / z2);

            assert!(
                f64::abs(eq1.re - 0_f64) < eps,
                "result: {}, exp: {}",
                eq1,
                0_f64
            );
            assert!(
                f64::abs(eq2.re + pn * (1_f64 + mu_0)) < eps,
                "result: {}, exp: {}",
                eq2,
                -pn * (1_f64 + mu_0)
            );
            assert!(
                f64::abs(eq3.re - 0_f64) < eps,
                "result: {}, exp: {}",
                eq3,
                0_f64
            );
            assert!(
                f64::abs(eq4.re - 0_f64) < eps,
                "result: {}, exp: {}",
                eq4,
                0_f64
            );
        }
    }
}
