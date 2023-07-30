use crate::{integration::definite_integral, utils::sum_calc};
use std::f64::consts::PI;

// returns c1, c2, c3, c4 coefficients
fn coefficients<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
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
    let coef = pn * 4_f64 * alpha * (1_f64 + mu_0);
    let e1 = f64::exp(alpha * b);
    let e2 = f64::exp(-alpha * b);

    let a1 = e1 * (2_f64 * b * alpha * alpha * mu_0 + 2_f64 * alpha * mu_0 + 2_f64 * alpha)
        - e2 * (-2_f64 * b * alpha * alpha * mu_0 + 2_f64 * alpha + 2_f64 * alpha * mu_0);
    let a2 = e1 * (2_f64 * b * alpha * alpha * mu_0 - 2_f64 * alpha)
        + e2 * (2_f64 * b * alpha * alpha * mu_0 + 2_f64 * alpha);
    let a3 = e1
        * (-2_f64 * g * b * alpha * alpha * mu_0 - 2_f64 * g * alpha * mu_0
            + 2_f64 * lambda * alpha)
        - e2 * (-2_f64 * g * b * alpha * alpha * mu_0 + 2_f64 * g * alpha * mu_0
            - 2_f64 * lambda * alpha);
    let a4 = e1 * (-2_f64 * g * b * alpha * alpha * mu_0 + (2_f64 * g + lambda) * 2_f64 * alpha)
        + e2 * (2_f64 * g * b * alpha * alpha * mu_0 + (2_f64 * g + lambda) * 2_f64 * alpha);

    let c1 = coef * a2 / (a4 * a1 - a2 * a3);
    let c2 = -coef * a1 / (a4 * a1 - a2 * a3);
    let c3 = -c1;
    let c4 = c2;

    (c1, c2, c3, c4)
}

fn function_un<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    y: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let coef = 1_f64 / alpha / 4_f64;
    let e1 = f64::exp(alpha * y);
    let e2 = f64::exp(-alpha * y);

    let (c1, c2, c3, c4) = coefficients(a, b, alpha, mu_0, g, lambda, load_function, eps);

    let res = c1 * e1 * (y * alpha * mu_0 + 2_f64 + mu_0)
        + c2 * e1 * (y * alpha * mu_0)
        + c3 * e2 * (y * alpha * mu_0 - 2_f64 - mu_0)
        + c4 * e2 * (-y * alpha * mu_0);

    res * coef
}

fn function_vn<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    y: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let coef = 1_f64 / alpha / 4_f64;
    let e1 = f64::exp(alpha * y);
    let e2 = f64::exp(-alpha * y);

    let (c1, c2, c3, c4) = coefficients(a, b, alpha, mu_0, g, lambda, load_function, eps);

    let res = c1 * e1 * (-y * alpha * mu_0)
        + c2 * e1 * (-y * alpha * mu_0 + 2_f64 + mu_0)
        + c3 * e2 * (y * alpha * mu_0)
        + c4 * e2 * (-y * alpha * mu_0 - 2_f64 - mu_0);

    res * coef
}

fn function_derivative_vn<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    y: f64,
    alpha: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let coef = 1.0 / alpha / 4.0;
    let e1 = f64::exp(alpha * y);
    let e2 = f64::exp(-alpha * y);

    let (c1, c2, c3, c4) = coefficients(a, b, alpha, mu_0, g, lambda, load_function, eps);

    let res = c1 * e1 * (-y * alpha * alpha * mu_0 - alpha * mu_0)
        + c2 * e1 * (-y * alpha * alpha * mu_0 + 2_f64 * alpha)
        + c3 * e2 * (-y * alpha * alpha * mu_0 + alpha * mu_0)
        + c4 * e2 * (y * alpha * alpha * mu_0 + 2_f64 * alpha);

    res * coef
}

fn function_u<F: Fn(f64) -> f64 + Send + Sync>(
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
    let n = 10;
    let start = 1;
    let f = |i| {
        if x != a && x != 0_f64 {
            let alpha = PI * i as f64 / a;
            2_f64
                * function_un(a, b, y, alpha, mu_0, g, lambda, load_function, eps)
                * f64::sin(alpha * x)
                / a
                / (1_f64 + mu_0)
        } else {
            0_f64
        }
    };

    sum_calc(&f, eps, start, n)
}

fn function_derivative_u_x<F: Fn(f64) -> f64 + Send + Sync>(
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
    let n = 10;
    let start = 1;
    let f = |i| {
        if x != a && x != 0_f64 {
            let alpha = PI * i as f64 / a;
            2_f64
                * alpha
                * function_un(a, b, y, alpha, mu_0, g, lambda, load_function, eps)
                * f64::cos(alpha * x)
                / a
                / (1_f64 + mu_0)
        } else {
            0_f64
        }
    };

    sum_calc(&f, eps, start, n)
}

fn function_v<F: Fn(f64) -> f64 + Send + Sync>(
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
    let p0 = definite_integral(0_f64, a, 100, eps, load_function);
    let v0 = -p0 * y / a / (2_f64 * g + lambda);
    let n = 10;
    let start = 1;
    let f = |i| {
        let alpha = PI * i as f64 / a;
        2_f64
            * function_vn(a, b, y, alpha, mu_0, g, lambda, load_function, eps)
            * f64::cos(alpha * x)
            / a
            / (1_f64 + mu_0)
    };

    v0 + sum_calc(&f, eps, start, n)
}

fn function_derivative_v_y<F: Fn(f64) -> f64 + Send + Sync>(
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
    let p0 = definite_integral(0_f64, a, 100, eps, load_function);
    let v0 = -p0 * y / a / (2_f64 * g + lambda);
    let start = 1;
    let n = 10;
    let f = |i| {
        let alpha = PI * i as f64 / a;
        2.0 * function_derivative_vn(a, b, y, alpha, mu_0, g, lambda, load_function, eps)
            * f64::cos(alpha * x)
            / a
            / (1_f64 + mu_0)
    };

    v0 + sum_calc(&f, eps, start, n)
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
    let d_ux = function_derivative_u_x(a, b, x, y, mu_0, g, lambda, load_function, eps);
    let d_vy = function_derivative_v_y(a, b, x, y, mu_0, g, lambda, load_function, eps);

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
    let d_ux = function_derivative_u_x(a, b, x, y, mu_0, g, lambda, load_function, eps);
    let d_vy = function_derivative_v_y(a, b, x, y, mu_0, g, lambda, load_function, eps);

    2_f64 * g * d_vy + lambda * d_vy + lambda * d_ux
}

#[cfg(test)]
mod run {
    use super::*;
    use crate::utils::{function_calculation, function_plot, g, lambda, mu_0, save_function};

    #[test]
    fn function_u_run() {
        let b = 15.0;
        // steel
        let puasson_coef = 0.25;
        let young_modulus = 200.0;

        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);

        let load_function = |x| x * x - 7.5;
        let eps = 0.1;

        let y = b;
        let a = b;

        let (x, y) = function_calculation(0.0, a, 100, |x| {
            function_sigma_y(a, b, x, y, mu_0, g, lambda, &load_function, eps)
        });

        let file_name = "static_1 function_u.txt";
        let file = save_function(&x, &y, "x", "u(x,y)", file_name);
        function_plot(&file)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::{g, lambda, mu_0};

    #[test]
    fn coefficients_test() {
        let load_function = |x| x * x;
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;
        let eps = 0.001;

        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);

        for i in 1..10 {
            let alpha = PI * i as f64 / a;

            let e1 = f64::exp(alpha * b);
            let e2 = f64::exp(-alpha * b);
            let pn = definite_integral(0_f64, a, 100, eps, &|x| {
                load_function(x) * f64::cos(alpha * x)
            });

            let (c1, c2, c3, c4) = coefficients(a, b, alpha, mu_0, g, lambda, &load_function, eps);

            let eq1 = c1
                * e1
                * (2_f64 * b * alpha * alpha * mu_0 + 2_f64 * alpha * mu_0 + 2_f64 * alpha)
                + c2 * e1 * (2_f64 * b * alpha * alpha * mu_0 - 2_f64 * alpha)
                + c3 * e2
                    * (-2_f64 * b * alpha * alpha * mu_0 + 2_f64 * alpha + 2_f64 * alpha * mu_0)
                + c4 * e2 * (2_f64 * b * alpha * alpha * mu_0 + 2_f64 * alpha);
            let eq2 = c1
                * e1
                * (-2_f64 * g * b * alpha * alpha * mu_0 - 2_f64 * g * alpha * mu_0
                    + 2_f64 * lambda * alpha)
                + c2 * e1
                    * (-2_f64 * g * b * alpha * alpha * mu_0
                        + (2_f64 * g + lambda) * 2_f64 * alpha)
                + c3 * e2
                    * (-2_f64 * g * b * alpha * alpha * mu_0 + 2_f64 * g * alpha * mu_0
                        - 2_f64 * lambda * alpha)
                + c4 * e2
                    * (2_f64 * g * b * alpha * alpha * mu_0 + (2_f64 * g + lambda) * 2_f64 * alpha);
            let eq3 = c1 * (2_f64 * alpha + 2_f64 * alpha * mu_0)
                + c2 * (-2_f64 * alpha)
                + c3 * (2_f64 * alpha + 2_f64 * alpha * mu_0)
                + c4 * (2_f64 * alpha);
            let eq4 = c2 * (2_f64 + mu_0) + c4 * (-2_f64 - mu_0);

            assert!(
                f64::abs(eq1 - 0_f64) < eps,
                "result: {}, exp: {}",
                eq1,
                0_f64
            );
            assert!(
                f64::abs(eq2 + pn * 4_f64 * alpha * (1_f64 + mu_0)) < eps,
                "result: {}, exp: {}",
                eq2,
                -pn * 4_f64 * alpha * (1_f64 + mu_0)
            );
            assert!(
                f64::abs(eq3 - 0_f64) < eps,
                "result: {}, exp: {}",
                eq3,
                0_f64
            );
            assert!(
                f64::abs(eq4 - 0_f64) < eps,
                "result: {}, exp: {}",
                eq4,
                0_f64
            );
        }
    }
}
