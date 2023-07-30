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

    let y = b;

    let a1 = e1
        * (alpha * alpha * mu_0 * y + 2.0 * alpha + mu_0 * alpha - alpha * (-alpha * mu_0 * y))
        - e2 * (-alpha * alpha * mu_0 * y + 2.0 * alpha + mu_0 * alpha
            - alpha * (alpha * mu_0 * y));

    let a2 = e1
        * (alpha * alpha * mu_0 * y + alpha * mu_0 - alpha * (-alpha * mu_0 * y + 2.0 + mu_0))
        + e2 * (alpha * alpha * mu_0 * y - alpha * mu_0 - alpha * (-alpha * mu_0 * y - 2.0 - mu_0));

    let a3 = e1
        * ((2.0 * g + lambda) * (-alpha * alpha * mu_0 * y - alpha * mu_0)
            + alpha * lambda * (alpha * mu_0 * y + 2.0 + mu_0))
        - e2 * ((2.0 * g + lambda) * (-alpha * alpha * mu_0 * y + alpha * mu_0)
            + alpha * lambda * (alpha * mu_0 * y - 2.0 - mu_0));

    let a4 = e1
        * ((2.0 * g + lambda) * (-alpha * alpha * mu_0 * y + 2.0 * alpha)
            + alpha * lambda * (alpha * mu_0 * y))
        + e2 * ((2.0 * g + lambda) * (alpha * alpha * mu_0 * y + 2.0 * alpha)
            + alpha * lambda * (-alpha * mu_0 * y));

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
    let coef = 1.0 / ((1.0 + mu_0) * 4.0 * alpha);
    let (c1, c2, c3, c4) = coefficients(a, b, alpha, mu_0, g, lambda, load_function, eps);

    let res = c1 * f64::exp(alpha * y) * (alpha * mu_0 * y + 2.0 + mu_0)
        + c2 * f64::exp(alpha * y) * (alpha * mu_0 * y)
        + c3 * f64::exp(-alpha * y) * (alpha * mu_0 * y - 2.0 - mu_0)
        + c4 * f64::exp(-alpha * y) * (-alpha * mu_0 * y);
    coef * res
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
    let coef = 1.0 / ((1.0 + mu_0) * 4.0 * alpha);

    let (c1, c2, c3, c4) = coefficients(a, b, alpha, mu_0, g, lambda, load_function, eps);

    let res = c1 * f64::exp(alpha * y) * (-alpha * mu_0 * y)
        + c2 * f64::exp(alpha * y) * (-alpha * mu_0 * y + 2.0 + mu_0)
        + c3 * f64::exp(-alpha * y) * (alpha * mu_0 * y)
        + c4 * f64::exp(-alpha * y) * (-alpha * mu_0 * y - 2.0 - mu_0);
    coef * res
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
    let coef = 1.0 / ((1.0 + mu_0) * 4.0 * alpha);

    let (c1, c2, c3, c4) = coefficients(a, b, alpha, mu_0, g, lambda, load_function, eps);

    let res = c1 * f64::exp(alpha * y) * (-alpha * alpha * mu_0 * y - alpha * mu_0)
        + c2 * f64::exp(alpha * y) * (-alpha * alpha * mu_0 * y + 2.0 * alpha)
        + c3 * f64::exp(-alpha * y) * (-alpha * alpha * mu_0 * y + alpha * mu_0)
        + c4 * f64::exp(-alpha * y) * (alpha * alpha * mu_0 * y + 2.0 * alpha);
    coef * res
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
        let alpha = PI * i as f64 / a;
        2_f64
            * function_un(a, b, y, alpha, mu_0, g, lambda, load_function, eps)
            * f64::sin(alpha * x)
            / a
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
        let alpha = PI * i as f64 / a;
        2_f64
            * alpha
            * function_un(a, b, y, alpha, mu_0, g, lambda, load_function, eps)
            * f64::cos(alpha * x)
            / a
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
    let v0 = -p0 / a / (2_f64 * g + lambda);
    let start = 1;
    let n = 10;
    let f = |i| {
        let alpha = PI * i as f64 / a;
        2.0 * function_derivative_vn(a, b, y, alpha, mu_0, g, lambda, load_function, eps)
            * f64::cos(alpha * x)
            / a
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

    (2.0 * g + lambda) * d_vy + lambda * d_ux
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

        let load_function = |x: f64| (x - b / 2.0) * (x - b / 2.0);
        let eps = 0.1;

        let y = b;
        let a = b;

        let (x, y) = function_calculation(0.0, a, 100, |x| {
            function_sigma_x(a, b, x, y, mu_0, g, lambda, &load_function, eps)
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
    use nalgebra::{Matrix2, Matrix2x1, RowVector1, RowVector2};

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

            let pn = definite_integral(0_f64, a, 100, eps, &|x| {
                load_function(x) * f64::cos(alpha * x)
            });
            println!("pn = {}", pn);
            let coef = 1.0 / ((1.0 + mu_0) * 4.0 * alpha);

            let (c1, c2, c3, c4) = coefficients(a, b, alpha, mu_0, g, lambda, &load_function, eps);

            let z1 = |y: f64| {
                let res = c1 * f64::exp(alpha * y) * (alpha * mu_0 * y + 2.0 + mu_0)
                    + c2 * f64::exp(alpha * y) * (alpha * mu_0 * y)
                    + c3 * f64::exp(-alpha * y) * (alpha * mu_0 * y - 2.0 - mu_0)
                    + c4 * f64::exp(-alpha * y) * (-alpha * mu_0 * y);
                coef * res
            };
            let z2 = |y: f64| {
                let res = c1 * f64::exp(alpha * y) * (-alpha * mu_0 * y)
                    + c2 * f64::exp(alpha * y) * (-alpha * mu_0 * y + 2.0 + mu_0)
                    + c3 * f64::exp(-alpha * y) * (alpha * mu_0 * y)
                    + c4 * f64::exp(-alpha * y) * (-alpha * mu_0 * y - 2.0 - mu_0);
                coef * res
            };

            let der_z1 = |y: f64| {
                let res = c1
                    * f64::exp(alpha * y)
                    * (alpha * alpha * mu_0 * y + 2.0 * alpha + mu_0 * alpha)
                    + c2 * f64::exp(alpha * y) * (alpha * alpha * mu_0 * y + alpha * mu_0)
                    + c3 * f64::exp(-alpha * y)
                        * (-alpha * alpha * mu_0 * y + 2.0 * alpha + mu_0 * alpha)
                    + c4 * f64::exp(-alpha * y) * (alpha * alpha * mu_0 * y - alpha * mu_0);
                coef * res
            };
            let der_z2 = |y: f64| {
                let res = c1 * f64::exp(alpha * y) * (-alpha * alpha * mu_0 * y - alpha * mu_0)
                    + c2 * f64::exp(alpha * y) * (-alpha * alpha * mu_0 * y + 2.0 * alpha)
                    + c3 * f64::exp(-alpha * y) * (-alpha * alpha * mu_0 * y + alpha * mu_0)
                    + c4 * f64::exp(-alpha * y) * (alpha * alpha * mu_0 * y + 2.0 * alpha);
                coef * res
            };

            let a_1 = Matrix2::from_rows(&[
                RowVector2::new(1.0, 0.0),
                RowVector2::new(0.0, 2.0 * g + lambda),
            ]);
            let a_2 = Matrix2::from_rows(&[
                RowVector2::new(0.0, -alpha),
                RowVector2::new(alpha * lambda, 0.0),
            ]);

            let res = a_1
                * Matrix2x1::from_rows(&[RowVector1::new(der_z1(b)), RowVector1::new(der_z2(b))])
                + a_2 * Matrix2x1::from_rows(&[RowVector1::new(z1(b)), RowVector1::new(z2(b))]);
            println!("U_0[Z_n]: {}", res);

            let a_1 = Matrix2::from_rows(&[RowVector2::new(1.0, 0.0), RowVector2::new(0.0, 0.0)]);
            let a_2 =
                Matrix2::from_rows(&[RowVector2::new(0.0, -alpha), RowVector2::new(0.0, 1.0)]);

            let res = a_1
                * Matrix2x1::from_rows(&[
                    RowVector1::new(der_z1(0.0)),
                    RowVector1::new(der_z2(0.0)),
                ])
                + a_2 * Matrix2x1::from_rows(&[RowVector1::new(z1(0.0)), RowVector1::new(z2(0.0))]);
            println!("U_1[Z_n]: {}", res);
        }
    }
}
