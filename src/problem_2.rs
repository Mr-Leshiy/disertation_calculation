use crate::{
    integration::{definite_integral, sqrt_gauss_integral, sqrt_gauss_integral_finit},
    matrices::{system_solve, Matrix},
    polynomials::chebyshev,
    utils::{
        function_calculation, g, lambda, mu_0, save_static, sum_calc_finit, surface_static_plot,
    },
    FunctionType, LoadFunction,
};
use clap::Parser;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use std::f64::consts::PI;

#[derive(Parser)]
#[clap(rename_all = "snake-case")]
pub struct Problem2 {
    #[clap(long)]
    a: f64,
    #[clap(long)]
    b: f64,
    #[clap(long)]
    n_x: u32,
    #[clap(long)]
    n_y: u32,
    #[clap(long)]
    puasson_coef: f64,
    #[clap(long)]
    young_modulus: f64,
    #[clap(long)]
    eps: f64,
    #[clap(long)]
    function_type: FunctionType,
    #[clap(long)]
    load_function: LoadFunction,
}

impl Problem2 {
    pub fn exec(self) {
        assert!(self.a < self.b);

        let mu_0 = mu_0(self.puasson_coef);
        let g = g(self.puasson_coef, self.young_modulus);
        let lambda = lambda(self.puasson_coef, self.young_modulus);

        let (x, y, z) = match self.function_type {
            FunctionType::U => {
                function_calculation(self.a, self.b, self.n_x, self.n_y, None, |x, y, _| {
                    function_u(
                        self.a,
                        self.b,
                        x,
                        y,
                        mu_0,
                        g,
                        lambda,
                        &|x| self.load_function.call(x),
                        self.eps,
                    )
                })
            }
            FunctionType::V => {
                function_calculation(self.a, self.b, self.n_x, self.n_y, None, |x, y, _| {
                    function_v(
                        self.a,
                        self.b,
                        x,
                        y,
                        mu_0,
                        g,
                        lambda,
                        &|x| self.load_function.call(x),
                        self.eps,
                    )
                })
            }
            FunctionType::SigmaX => {
                function_calculation(self.a, self.b, self.n_x, self.n_y, None, |x, y, _| {
                    function_sigma_x(
                        self.a,
                        self.b,
                        x,
                        y,
                        mu_0,
                        g,
                        lambda,
                        &|x| self.load_function.call(x),
                        self.eps,
                    )
                })
            }
            FunctionType::SigmaY => {
                function_calculation(self.a, self.b, self.n_x, self.n_y, None, |x, y, _| {
                    function_sigma_y(
                        self.a,
                        self.b,
                        x,
                        y,
                        mu_0,
                        g,
                        lambda,
                        &|x| self.load_function.call(x),
                        self.eps,
                    )
                })
            }
        };

        let file_name = format!("problem_2, a: {0}, b: {1}, n_x: {2}, n_y: {3}, puasson_coef: {4}, young_modulus: {5}, eps: {6}, function_type: {7}, load_function: {8}", self.a, self.b, self.n_x, self.n_y, self.puasson_coef, self.young_modulus, self.eps, self.function_type, self.load_function);
        let path = save_static(&x, &y, &z[0], &file_name);
        surface_static_plot(&path);
    }
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
    let coef = 1.0 + mu_0;
    let e_1 = f64::exp(-alpha * b);
    let e_2 = e_1 * e_1;
    let inv_alpha = 1.0 / alpha;

    let x1 = b * mu_0 + inv_alpha + inv_alpha * mu_0;
    let x2 = b * mu_0 - inv_alpha;
    let x3 = -b * mu_0 + inv_alpha + inv_alpha * mu_0;
    let x4 = b * mu_0 + inv_alpha;
    let x5 = -g * b * mu_0 - inv_alpha * g * mu_0 + inv_alpha * lambda;
    let x6 = -g * b * mu_0 + 2.0 * inv_alpha * g + inv_alpha * lambda;
    let x7 = -g * b * mu_0 + inv_alpha * g * mu_0 - inv_alpha * lambda;
    let x8 = g * b * mu_0 + inv_alpha * 2.0 * g + inv_alpha * lambda;
    let x9 = 1.0 + mu_0;
    let x10 = -1.0;
    let x11 = 2.0 + mu_0;

    let y1 = x2 + x4 * e_2;
    let y2 = -x1 + x3 * e_2;
    let y3 = x6 + x8 * e_2;
    let y4 = -x5 + x7 * e_2;
    let y5 = -(g + lambda) * x11
        + alpha * alpha * e_2 * (x3 * x6 - x1 * x8 - x2 * x7 + x4 * x5)
        + alpha * alpha * e_2 * e_2 * (x3 * x8 - x4 * x7);

    let d_1_0 = -coef * 2.0 * y3 / y5;
    let d_2_0 = coef * 2.0 * y1 / y5;
    let d_3_0 = -coef * 2.0 * y4 / y5;
    let d_4_0 = coef * 2.0 * y2 / y5;

    let f_1_0 = coef * 2.0 * y3 / y5;
    let f_2_0 = -coef * 2.0 * y1 / y5;
    let f_3_0 = -coef * 2.0 * y4 / y5;
    let f_4_0 = coef * 2.0 * y2 / y5;

    let z1 = (g + lambda) * x11 + alpha * alpha * e_2 * (x1 * x8 - x4 * x5);
    let z2 = x3 * x6 - x2 * x7 + e_2 * e_2 * (x3 * x8 - x4 * x7);
    let z3 = x1 * x7 - x3 * x5;
    let z4 =
        (g + lambda) * x11 + e_2 * (x4 * x5 * x10 - x4 * x6 * x9 - x1 * x8 * x10 + x2 * x8 * x9);
    let z5 = x2 * x8 * x9 - x4 * x6 * x9 + x2 * x7 * x10 - x3 * x6 * x10
        + alpha * alpha * e_2 * e_2 * (x4 * x7 * x10 - x3 * x8 * x10);
    let z6 = (g + lambda) * x11 * coef
        + e_2 * (x3 * x5 * x10 - x3 * x6 * x9 + x2 * x7 * x9 - x1 * x7 * x10);
    let z7 = x3 * x5 * x10 - x1 * x7 * x10 - x1 * x8 * x9
        + x4 * x5 * x9
        + alpha * alpha * e_2 * e_2 * (x3 * x8 * x9 - x4 * x7 * x9);

    let d_1_1 = coef * 2.0 * z2 / (x9 * y5);
    let d_2_1 = coef * 4.0 * z5 / (x9 * x11 * y5);
    let d_3_1 = coef * 2.0 * z3 / (x9 * y5);
    let d_4_1 = coef * 4.0 * z7 / (x9 * x11 * y5);

    let f_1_1 = -coef * 2.0 * z1 / (x9 * y5);
    let f_2_1 = -coef * 4.0 * z4 / (x9 * x11 * y5);
    let f_3_1 = coef * 2.0 * z3 / (x9 * y5);
    let f_4_1 = coef * 4.0 * z6 / (x9 * x11 * y5);

    (
        (d_1_0, d_2_0, d_3_0, d_4_0, f_1_0, f_2_0, f_3_0, f_4_0),
        (d_1_1, d_2_1, d_3_1, d_4_1, f_1_1, f_2_1, f_3_1, f_4_1),
    )
}

#[allow(dead_code)]
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

    let d_1_0 = e * alpha * d_1_0;
    let d_2_0 = e * alpha * d_2_0;
    let d_3_0 = e * alpha * d_3_0;
    let d_4_0 = e * alpha * d_4_0;

    let f_1_0 = e * alpha * f_1_0;
    let f_2_0 = e * alpha * f_2_0;
    let f_3_0 = e * alpha * f_3_0;
    let f_4_0 = e * alpha * f_4_0;

    let d_1_1 = alpha * alpha * e * e * d_1_1;
    let d_2_1 = alpha * alpha * alpha * e * e * d_2_1;
    let d_3_1 = alpha * alpha * e * e * d_3_1;
    let d_4_1 = alpha * alpha * alpha * e * e * d_4_1;

    let f_1_1 = f_1_1;
    let f_2_1 = alpha * f_2_1;
    let f_3_1 = alpha * alpha * e * e * f_3_1;
    let f_4_1 = alpha * f_4_1;

    (
        (d_1_0, d_2_0, d_3_0, d_4_0, f_1_0, f_2_0, f_3_0, f_4_0),
        (d_1_1, d_2_1, d_3_1, d_4_1, f_1_1, f_2_1, f_3_1, f_4_1),
    )
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
    let e1: f64 = f64::exp(-2.0 * alpha * x);
    let e2 = f64::exp(2.0 * alpha * x);
    let e3 = f64::exp(-2.0 * alpha * b);

    let (
        (d_1_0, d_2_0, d_3_0, d_4_0, f_1_0, f_2_0, f_3_0, f_4_0),
        (d_1_1, d_2_1, d_3_1, d_4_1, f_1_1, f_2_1, f_3_1, f_4_1),
    ) = coefficients_1(b, alpha, mu_0, g, lambda);

    let psi_1_0 = ((x * alpha * mu_0 + 2.0 + mu_0) * d_1_0 + x * alpha * mu_0 * d_3_0)
        + e1 * ((x * alpha * mu_0 - 2.0 - mu_0) * f_1_0 - x * alpha * mu_0 * f_3_0);
    let psi_2_0 = ((x * alpha * mu_0 + 2.0 + mu_0) * d_2_0 + x * alpha * mu_0 * d_4_0)
        + e1 * ((x * alpha * mu_0 - 2.0 - mu_0) * f_2_0 - x * alpha * mu_0 * f_4_0);
    let psi_3_0 = (-x * alpha * mu_0 * d_1_0 + (-x * alpha * mu_0 + 2.0 + mu_0) * d_3_0)
        + e1 * (x * alpha * mu_0 * f_1_0 + (-x * alpha * mu_0 - 2.0 - mu_0) * f_3_0);
    let psi_4_0 = (-x * alpha * mu_0 * d_2_0 + (-x * alpha * mu_0 + 2.0 + mu_0) * d_4_0)
        + e1 * (x * alpha * mu_0 * f_2_0 + (-x * alpha * mu_0 - 2.0 - mu_0) * f_4_0);

    let psi_1_1 = (x * mu_0 - 2.0 / alpha - mu_0 / alpha) * f_1_1
        + e2 * e3
            * ((x * alpha * mu_0 + 2.0 + mu_0) * alpha * d_1_1 + x * alpha * mu_0 * alpha * d_3_1)
        + e3 * (-x * alpha * mu_0 * alpha * f_3_1);
    let psi_2_1 = (x * mu_0 - 2.0 / alpha - mu_0 / alpha) * f_2_1 - x * mu_0 * f_4_1
        + e2 * e3
            * ((x * alpha * mu_0 + 2.0 + mu_0) * alpha * d_2_1 + x * alpha * alpha * mu_0 * d_4_1);
    let psi_3_1 = -x * mu_0 * f_1_1
        + e2 * e3
            * (-x * alpha * alpha * mu_0 * d_1_1
                + (-x * alpha * mu_0 + 2.0 + mu_0) * alpha * d_3_1)
        + e3 * (-x * alpha * mu_0 - 2.0 - mu_0) * alpha * f_3_1;
    let psi_4_1 = x * mu_0 * f_2_1
        + (-x * mu_0 - 2.0 / alpha - mu_0 / alpha) * f_4_1
        + e2 * e3
            * (-x * alpha * alpha * mu_0 * d_2_1
                + (-x * alpha * mu_0 + 2.0 + mu_0) * alpha * d_4_1);

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
    let e1: f64 = f64::exp(-2.0 * alpha * x);
    let e2 = f64::exp(2.0 * alpha * x);
    let e3 = f64::exp(-2.0 * alpha * b);

    let (
        (d_1_0, d_2_0, d_3_0, d_4_0, f_1_0, f_2_0, f_3_0, f_4_0),
        (d_1_1, d_2_1, d_3_1, d_4_1, f_1_1, f_2_1, f_3_1, f_4_1),
    ) = coefficients_1(b, alpha, mu_0, g, lambda);

    let psi_1_0 = ((x * alpha * mu_0 + 2.0 + 2.0 * mu_0) * d_1_0
        + (x * alpha * mu_0 + mu_0) * d_3_0)
        + e1 * ((-x * alpha * mu_0 + 2.0 + 2.0 * mu_0) * f_1_0 + (x * alpha * mu_0 - mu_0) * f_3_0);
    let psi_2_0 = ((x * alpha * mu_0 + 2.0 + 2.0 * mu_0) * d_2_0
        + (x * alpha * mu_0 + mu_0) * d_4_0)
        + e1 * ((-x * alpha * mu_0 + 2.0 + 2.0 * mu_0) * f_2_0 + (x * alpha * mu_0 - mu_0) * f_4_0);
    let psi_3_0 = ((-x * alpha * mu_0 - mu_0) * d_1_0 + (-x * alpha * mu_0 + 2.0) * d_3_0)
        + e1 * ((-x * alpha * mu_0 + mu_0) * f_1_0 + (x * alpha * mu_0 + 2.0) * f_3_0);
    let psi_4_0 = ((-x * alpha * mu_0 - mu_0) * d_2_0 + (-x * alpha * mu_0 + 2.0) * d_4_0)
        + e1 * ((-x * alpha * mu_0 + mu_0) * f_2_0 + (x * alpha * mu_0 + 2.0) * f_4_0);

    let psi_1_1 = (-x * mu_0 + 2.0 / alpha + 2.0 * mu_0 / alpha) * f_1_1
        + e2 * e3
            * ((x * alpha * mu_0 + 2.0 + 2.0 * mu_0) * alpha * d_1_1
                + (x * alpha * mu_0 + mu_0) * alpha * d_3_1)
        + e3 * ((x * alpha * mu_0 - mu_0) * f_3_1);
    let psi_2_1 = (-x * mu_0 + 2.0 / alpha + 2.0 * mu_0 / alpha) * f_2_1
        + (x * mu_0 - mu_0 / alpha) * f_4_1
        + e2 * e3
            * ((x * alpha * mu_0 + 2.0 + 2.0 * mu_0) * alpha * d_2_1
                + (x * alpha * mu_0 + mu_0) * alpha * d_4_1);
    let psi_3_1 = (-x * mu_0 + mu_0 / alpha) * f_1_1
        + e2 * e3
            * ((-x * alpha * mu_0 - mu_0) * d_1_1 + (-x * alpha * mu_0 + 2.0) * alpha * d_3_1)
        + e3 * ((x * alpha * mu_0 + 2.0) * alpha * f_3_1);
    let psi_4_1 = (-x * mu_0 + mu_0 / alpha) * f_2_1
        + (x * mu_0 + 2.0 / alpha) * f_4_1
        + e2 * e3
            * ((-x * alpha * mu_0 - mu_0) * alpha * d_2_1
                + (-x * alpha * mu_0 + 2.0) * alpha * d_4_1);

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
    let e1 = f64::exp(-alpha * x);
    let e2 = f64::exp(alpha * (x - b));

    let ((psi_1_0, psi_2_0, psi_3_0, psi_4_0), (psi_1_1, psi_2_1, psi_3_1, psi_4_1)) =
        psi_1(x, b, alpha, mu_0, g, lambda);

    let psi_1_0 = e2 * psi_1_0;
    let psi_2_0 = e2 * psi_2_0;
    let psi_3_0 = e2 * psi_3_0;
    let psi_4_0 = e2 * psi_4_0;

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
    let e1 = f64::exp(-alpha * x);
    let e2 = f64::exp(alpha * (x - b));

    let ((psi_1_0, psi_2_0, psi_3_0, psi_4_0), (psi_1_1, psi_2_1, psi_3_1, psi_4_1)) =
        der_psi_1(x, b, alpha, mu_0, g, lambda);

    let psi_1_0 = alpha * e2 * psi_1_0;
    let psi_2_0 = alpha * e2 * psi_2_0;
    let psi_3_0 = alpha * e2 * psi_3_0;
    let psi_4_0 = alpha * e2 * psi_4_0;

    let psi_1_1 = alpha * e1 * psi_1_1;
    let psi_2_1 = alpha * alpha * e1 * psi_2_1;
    let psi_3_1 = alpha * e1 * psi_3_1;
    let psi_4_1 = alpha * alpha * e1 * psi_4_1;

    (
        (psi_1_0, psi_2_0, psi_3_0, psi_4_0),
        (psi_1_1, psi_2_1, psi_3_1, psi_4_1),
    )
}

fn a_1<F: Fn(f64) -> f64 + Send + Sync>(a: f64, x: f64, load_function: &F, eps: f64) -> f64 {
    const N: i32 = 10;
    let p0 = definite_integral(0_f64, a, 30, eps, load_function);
    let sum_fn = |t: f64| {
        if t == a {
            return 0.0;
        }
        let sum = (1..N)
            .into_par_iter()
            .map(|n| {
                let alpha = PI / a * n as f64;
                f64::cos(alpha * t) * f64::cos(alpha * x)
            })
            .sum::<f64>();
        sum
    };
    let sum = if x != a {
        definite_integral(0.0, a, 10, eps, &sum_fn)
    } else {
        0.0
    };

    let res = 2.0 / a * sum - load_function(x) + p0 / a;
    res
}

fn a_2(a: f64, b: f64, mu_0: f64, g: f64, lambda: f64) -> f64 {
    2.0 * a * g * b * b * mu_0 * mu_0 / (PI * (1.0 + mu_0) * (g + lambda))
}

fn a_3(a: f64, b: f64, x: f64, ksi: f64, mu_0: f64, g: f64, lambda: f64) -> f64 {
    const N: i32 = 10;
    let sum1 = (1..N)
        .into_par_iter()
        .map(|n| {
            let alpha = PI / a * n as f64;
            let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
            let a1 = sign * f64::exp(alpha * (ksi - 2.0 * b)) * f64::cos(alpha * x) / alpha;
            let ((_, psi_2_0_ksi, _, psi_4_0_ksi), _) = psi_1(ksi, b, alpha, mu_0, g, lambda);
            let (_, (psi_1_1_b, psi_2_1_b, _, _)) = psi_1(b, b, alpha, mu_0, g, lambda);
            let (_, (_, _, der_psi_3_1_b, der_psi_4_1_b)) = der_psi_1(b, b, alpha, mu_0, g, lambda);
            let a2 = (der_psi_3_1_b * psi_2_0_ksi / alpha + der_psi_4_1_b * psi_4_0_ksi)
                * (2.0 * g + lambda)
                + (psi_1_1_b * psi_2_0_ksi / alpha + psi_2_1_b * psi_4_0_ksi) * lambda;

            a1 * a2
        })
        .sum::<f64>();
    let sum2 = (0..N)
        .into_par_iter()
        .map(|n| {
            let alpha = 2.0 * n as f64 + 1.0;
            let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
            let a1 = sign
                * f64::exp(alpha * (ksi - 2.0 * b) * PI / (2.0 * a))
                * f64::cos(alpha * x * PI / (2.0 * a))
                / alpha;
            a1
        })
        .sum::<f64>();
    sum1 - a_2(a, b, mu_0, g, lambda) * sum2
}

fn a_4(a: f64, b: f64, h: f64) -> f64 {
    f64::sinh(f64::acosh(
        h * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a))) + f64::cosh(PI * b / (2.0 * a)),
    ))
}

fn a_5(a: f64, b: f64, mu_0: f64) -> f64 {
    -2.0 * (1.0 + mu_0) * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
}

fn f_m<F: Fn(f64) -> f64 + Send + Sync>(
    m: usize,
    a: f64,
    b: f64,
    mu_0: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let f = |l| {
        let x = 2.0 * a / PI * f64::acos(l);
        let a_1 = a_1(a, x, load_function, eps);
        let a_5 = a_5(a, b, mu_0);
        let cheb = chebyshev(l, 2 * m + 1);
        cheb * a_1 / a_5
    };
    sqrt_gauss_integral(10, eps, &f)
}

fn g_k_m(m: usize, k: usize, a: f64, b: f64, mu_0: f64, g: f64, lambda: f64, eps: f64) -> f64 {
    let f = |l| {
        let cheb = chebyshev(l, 2 * m + 1);
        let f = |h| {
            let h = f64::abs(h);
            let cheb = chebyshev(h, 2 * k + 1);
            let h = 2.0 * b
                - 2.0 * a / PI
                    * f64::acosh(
                        h * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
                            + f64::cosh(PI * b / (2.0 * a)),
                    );
            let l = 2.0 * a / PI * f64::acos(l);
            let a_3 = a_3(a, b, h, l, mu_0, g, lambda);
            let a_3 = 2.0 * a_3 + 2.0 * g + lambda;
            let a_2 = a_2(a, b, mu_0, g, lambda);
            a_3 * cheb / a_2
        };
        cheb * sqrt_gauss_integral(10, eps, &f) / 2.0
    };
    let res = sqrt_gauss_integral(10, eps, &f) / PI;
    res
}

fn phi<F: Fn(f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> Matrix {
    const N: usize = 3;
    let left: Vec<Vec<_>> = vec![(0..N)
        .into_par_iter()
        .map(|m| f_m(m as usize, a, b, mu_0, load_function, eps))
        .collect()];

    let g_k_m: Vec<Vec<_>> = (0..N)
        .into_par_iter()
        .map(|k| {
            (0..N)
                .into_par_iter()
                .map(|m| {
                    let res = g_k_m(m, k, a, b, mu_0, g, lambda, eps);
                    res
                })
                .collect()
        })
        .collect();

    let a: Vec<Vec<_>> = (0..N)
        .into_par_iter()
        .map(|m| {
            (0..N)
                .into_par_iter()
                .map(|k| {
                    if m == k {
                        g_k_m[k][m] + PI / 4.0 / (2.0 * m as f64 + 1_f64)
                    } else {
                        g_k_m[k][m]
                    }
                })
                .collect()
        })
        .collect();

    system_solve(Matrix::new(a), Matrix::new(left), eps)
}

fn unknown_function<F: Fn(f64) -> f64 + Send + Sync>(
    h: f64,
    a: f64,
    b: f64,
    mu_0: f64,
    g: f64,
    lambda: f64,
    load_function: &F,
    eps: f64,
) -> f64 {
    let phi = phi(a, b, mu_0, g, lambda, load_function, eps);

    let a2 = a_2(a, b, mu_0, g, lambda);
    let a4 = a_4(a, b, h);
    // let sqrt = f64::sqrt(1.0 - h * h);
    let f = |i| {
        if i < phi.m() {
            phi.get_element(0, i) * chebyshev(h, 2 * i + 1)
        } else {
            0.0
        }
    };
    let sum = sum_calc_finit(&f, 0, 10);

    a4 * sum / a2
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
    let pn = definite_integral(0_f64, a, 100, eps, &|x| {
        load_function(x) * f64::cos(alpha * x)
    });
    let ((y_psi_1_0, y_psi_2_0, _, __), (y_psi_1_1, y_psi_2_1, _, _)) =
        psi_2(y, b, alpha, mu_0, g, lambda);
    let f = |ksi: f64| {
        let ksi = f64::abs(ksi);
        let f_val = unknown_function(ksi, a, b, mu_0, g, lambda, load_function, eps);
        let a1 = f64::sqrt(
            ksi * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
                + f64::cosh(PI * b / (2.0 * a))
                - 1.0,
        ) * f64::sqrt(
            ksi * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
                + f64::cosh(PI * b / (2.0 * a))
                + 1.0,
        );
        let h = 2.0 * b
            - 2.0 * a / PI
                * f64::acosh(
                    ksi * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
                        + f64::cosh(PI * b / (2.0 * a)),
                );
        let ((_, h_psi_2_0, _, h_psi_4_0), (_, h_psi_2_1, _, h_psi_4_1)) =
            psi_2(h, b, alpha, mu_0, g, lambda);
        let g1 = if y < h {
            (y_psi_1_0 * h_psi_2_1 + y_psi_2_0 * h_psi_4_1) / alpha / alpha / alpha
        } else {
            (y_psi_1_1 * h_psi_2_0 + y_psi_2_1 * h_psi_4_0) / alpha / alpha / alpha
        };
        f_val * g1 / a1
    };
    let int_val = sqrt_gauss_integral_finit(4, &f) / 2.0;
    let coef =
        f64::cos(alpha * a) * 2.0 * a * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
            / PI;
    coef * int_val - y_psi_2_0 * pn
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
    let pn = definite_integral(0_f64, a, 100, eps, &|x| {
        load_function(x) * f64::cos(alpha * x)
    });
    let ((_, _, y_psi_3_0, y_psi_4_0), (_, _, y_psi_3_1, y_psi_4_1)) =
        psi_2(y, b, alpha, mu_0, g, lambda);
    let f = |ksi: f64| {
        let ksi = f64::abs(ksi);
        let f_val = unknown_function(ksi, a, b, mu_0, g, lambda, load_function, eps);
        let a1 = f64::sqrt(
            ksi * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
                + f64::cosh(PI * b / (2.0 * a))
                - 1.0,
        ) * f64::sqrt(
            ksi * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
                + f64::cosh(PI * b / (2.0 * a))
                + 1.0,
        );
        let h = 2.0 * b
            - 2.0 * a / PI
                * f64::acosh(
                    ksi * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
                        + f64::cosh(PI * b / (2.0 * a)),
                );
        let ((_, h_psi_2_0, _, h_psi_4_0), (_, h_psi_2_1, _, h_psi_4_1)) =
            psi_2(h, b, alpha, mu_0, g, lambda);
        let g3 = if y < h {
            (y_psi_3_0 * h_psi_2_1 + y_psi_4_0 * h_psi_4_1) / alpha / alpha / alpha
        } else {
            (y_psi_3_1 * h_psi_2_0 + y_psi_4_1 * h_psi_4_0) / alpha / alpha / alpha
        };

        g3 * f_val / a1
    };
    let int_val = sqrt_gauss_integral(20, eps, &f) / 2.0;
    let coef =
        f64::cos(alpha * a) * 2.0 * a * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
            / PI;

    coef * int_val * int_val - y_psi_4_0 * pn
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
    let pn = definite_integral(0.0, a, 100, eps, &|x| {
        load_function(x) * f64::cos(alpha * x)
    });
    let ((_, _, y_psi_3_0, y_psi_4_0), (_, _, y_psi_3_1, y_psi_4_1)) =
        der_psi_2(y, b, alpha, mu_0, g, lambda);
    let f = |ksi: f64| {
        let ksi = f64::abs(ksi);
        let f_val = unknown_function(ksi, a, b, mu_0, g, lambda, load_function, eps);
        let a1 = f64::sqrt(
            ksi * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
                + f64::cosh(PI * b / (2.0 * a))
                - 1.0,
        ) * f64::sqrt(
            ksi * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
                + f64::cosh(PI * b / (2.0 * a))
                + 1.0,
        );
        let h = 2.0 * b
            - 2.0 * a / PI
                * f64::acosh(
                    ksi * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
                        + f64::cosh(PI * b / (2.0 * a)),
                );
        let ((_, h_psi_2_0, _, h_psi_4_0), (_, h_psi_2_1, _, h_psi_4_1)) =
            psi_2(h, b, alpha, mu_0, g, lambda);
        let g3 = if y < h {
            (y_psi_3_0 * h_psi_2_1 + y_psi_4_0 * h_psi_4_1) / alpha / alpha / alpha
        } else {
            (y_psi_3_1 * h_psi_2_0 + y_psi_4_1 * h_psi_4_0) / alpha / alpha / alpha
        };

        g3 * f_val / a1
    };
    let int_val = sqrt_gauss_integral(20, eps, &f) / 2.0;
    let coef =
        f64::cos(alpha * a) * 2.0 * a * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
            / PI;

    coef * int_val * int_val - y_psi_4_0 * pn
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
    let n = 3;
    let start = 1;
    let f = |i| {
        if x != a && x != 0.0 {
            let alpha = PI / a * i as f64;
            2.0 * function_un(a, b, y, alpha, mu_0, g, lambda, load_function, eps)
                * f64::sin(alpha * x)
                / a
        } else {
            0.0
        }
    };

    sum_calc_finit(&f, start, n)
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
    let n = 3;
    let start = 1;
    let f = |i| {
        if x != a && x != 0_f64 {
            let alpha = PI / a * i as f64;
            2.0 * alpha
                * function_un(a, b, y, alpha, mu_0, g, lambda, load_function, eps)
                * f64::cos(alpha * x)
                / a
        } else {
            0_f64
        }
    };

    sum_calc_finit(&f, start, n)
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
    let p0 = definite_integral(0_f64, a, 30, eps, load_function);
    let f = |ksi: f64| {
        let ksi = f64::abs(ksi);
        let f_val = unknown_function(ksi, a, b, mu_0, g, lambda, load_function, eps);
        let a1 = f64::sqrt(
            ksi * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
                + f64::cosh(PI * b / (2.0 * a))
                - 1.0,
        ) * f64::sqrt(
            ksi * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
                + f64::cosh(PI * b / (2.0 * a))
                + 1.0,
        );
        let h = 2.0 * b
            - 2.0 * a / PI
                * f64::acosh(
                    ksi * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
                        + f64::cosh(PI * b / (2.0 * a)),
                );
        let g = if y < h { -y } else { -h };
        f_val * g / a1
    };
    let coef = 2.0 * a * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a))) / PI;
    let int_val = sqrt_gauss_integral_finit(4, &f) / 2.0;
    let v_0 = 2.0 * coef * int_val / a - y * p0 / a / (2.0 * g + lambda);

    let n = 3;
    let start = 1;
    let f = |i| {
        let alpha = PI / a * i as f64;
        2.0 * function_vn(a, b, y, alpha, mu_0, g, lambda, load_function, eps) * f64::cos(alpha * x)
            / a
    };

    v_0 + sum_calc_finit(&f, start, n)
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
    let p0 = definite_integral(0_f64, a, 30, eps, load_function);
    let f = |ksi: f64| {
        let ksi = f64::abs(ksi);
        let f_val = unknown_function(ksi, a, b, mu_0, g, lambda, load_function, eps);
        let a1 = f64::sqrt(
            ksi * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
                + f64::cosh(PI * b / (2.0 * a))
                - 1.0,
        ) * f64::sqrt(
            ksi * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
                + f64::cosh(PI * b / (2.0 * a))
                + 1.0,
        );
        let h = 2.0 * b
            - 2.0 * a / PI
                * f64::acosh(
                    ksi * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a)))
                        + f64::cosh(PI * b / (2.0 * a)),
                );
        let g = if y < h { -1.0 } else { 0.0 };
        f_val * g / a1
    };
    let coef = 2.0 * a * (f64::cosh(PI * b / a) - f64::cosh(PI * b / (2.0 * a))) / PI;
    let int_val = sqrt_gauss_integral_finit(4, &f) / 2.0;
    let v_0 = 2.0 * coef * int_val / a - p0 / a / (2.0 * g + lambda);

    let n = 3;
    let start = 1;
    let f = |i| {
        let alpha = PI / a * i as f64;
        2.0 * function_derivative_vn(a, b, y, alpha, mu_0, g, lambda, load_function, eps)
            * f64::cos(alpha * x)
            / a
    };

    v_0 + sum_calc_finit(&f, start, n)
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

    2.0 * g * d_ux + lambda * d_vy + lambda * d_ux
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

    2.0 * g * d_vy + lambda * d_vy + lambda * d_ux
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
        let e1 = f64::exp(-2.0 * alpha * x);
        let e2 = f64::exp(2.0 * alpha * x);
        let e3 = f64::exp(-2.0 * alpha * b);

        let (
            (d_1_0, d_2_0, d_3_0, d_4_0, f_1_0, f_2_0, f_3_0, f_4_0),
            (d_1_1, d_2_1, d_3_1, d_4_1, f_1_1, f_2_1, f_3_1, f_4_1),
        ) = coefficients_1(b, alpha, mu_0, g, lambda);

        let psi_1_0 = ((x * alpha * mu_0 + 2.0 + 3.0 * mu_0) * d_1_0
            + (x * alpha * mu_0 + 2.0 * mu_0) * d_3_0)
            + e1 * ((x * alpha * mu_0 - 2.0 - 3.0 * mu_0) * f_1_0
                + (-x * alpha * mu_0 + 2.0 * mu_0) * f_3_0);
        let psi_2_0 = ((x * alpha * mu_0 + 2.0 + 3.0 * mu_0) * d_2_0
            + (x * alpha * mu_0 + 2.0 * mu_0) * d_4_0)
            + e1 * ((x * alpha * mu_0 - 2.0 - 3.0 * mu_0) * f_2_0
                + (-x * alpha * mu_0 + 2.0 * mu_0) * f_4_0);
        let psi_3_0 = ((-x * alpha * mu_0 - 2.0 * mu_0) * d_1_0
            + (-x * alpha * mu_0 + 2.0 - mu_0) * d_3_0)
            + e1 * ((x * alpha * mu_0 - 2.0 * mu_0) * f_1_0
                + (-x * alpha * mu_0 - 2.0 + mu_0) * f_3_0);
        let psi_4_0 = ((-x * alpha * mu_0 - 2.0 * mu_0) * d_2_0
            + (-x * alpha * mu_0 + 2.0 - mu_0) * d_4_0)
            + e1 * ((x * alpha * mu_0 - 2.0 * mu_0) * f_2_0
                + (-x * alpha * mu_0 - 2.0 + mu_0) * f_4_0);

        let psi_1_1 = (x * mu_0 - 2.0 / alpha - 3.0 * mu_0 / alpha) * f_1_1
            + e2 * e3
                * ((x * alpha * mu_0 + 2.0 + 3.0 * mu_0) * alpha * d_1_1
                    + (x * alpha * mu_0 + 2.0 * mu_0) * alpha * d_3_1)
            + e3 * ((-x * alpha * mu_0 + 2.0 * mu_0) * alpha * f_3_1);
        let psi_2_1 = (x * mu_0 - 2.0 / alpha - 3.0 * mu_0 / alpha) * f_2_1
            + (-x * mu_0 + 2.0 * mu_0 / alpha) * f_4_1
            + e2 * e3
                * ((x * alpha * mu_0 + 2.0 + 3.0 * mu_0) * alpha * d_2_1
                    + (x * alpha * mu_0 + 2.0 * mu_0) * alpha * d_4_1);
        let psi_3_1 = (x * mu_0 - 2.0 * mu_0 / alpha) * f_1_1
            + e2 * e3
                * ((-x * alpha * mu_0 - 2.0 * mu_0) * alpha * d_1_1
                    + (-x * alpha * mu_0 + 2.0 - mu_0) * alpha * d_3_1)
            + e3 * ((-x * alpha * mu_0 - 2.0 + mu_0) * alpha * f_3_1);
        let psi_4_1 = (x * mu_0 - 2.0 * mu_0 / alpha) * f_2_1
            + (-x * mu_0 - 2.0 / alpha + mu_0 / alpha) * f_4_1
            + e2 * e3
                * ((-x * alpha * mu_0 - 2.0 * mu_0) * alpha * d_2_1
                    + (-x * alpha * mu_0 + 2.0 - mu_0) * alpha * d_4_1);

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
        let e1 = f64::exp(-alpha * x);
        let e2 = f64::exp(alpha * (x - b));

        let ((psi_1_0, psi_2_0, psi_3_0, psi_4_0), (psi_1_1, psi_2_1, psi_3_1, psi_4_1)) =
            der_2_psi_1(x, b, alpha, mu_0, g, lambda);

        let psi_1_0 = alpha * alpha * e2 * psi_1_0;
        let psi_2_0 = alpha * alpha * e2 * psi_2_0;
        let psi_3_0 = alpha * alpha * e2 * psi_3_0;
        let psi_4_0 = alpha * alpha * e2 * psi_4_0;

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

        for i in 2..5 {
            println!("-----");
            let alpha = PI / a * i as f64;

            let e1 = f64::exp(alpha * b);
            let e2 = f64::exp(-alpha * b);

            let x1 = 2.0 * b * alpha * alpha * mu_0 + 2.0 * alpha + 2.0 * alpha * mu_0;
            let x2 = 2.0 * b * alpha * alpha * mu_0 - 2.0 * alpha;
            let x3 = -2.0 * b * alpha * alpha * mu_0 + 2.0 * alpha + 2.0 * alpha * mu_0;
            let x4 = 2.0 * b * alpha * alpha * mu_0 + 2.0 * alpha;
            let x5 =
                -2.0 * g * b * alpha * alpha * mu_0 - 2.0 * g * alpha * mu_0 + 2.0 * alpha * lambda;
            let x6 = -2.0 * g * b * alpha * alpha * mu_0 + 4_f64 * g * alpha + 2.0 * alpha * lambda;
            let x7 =
                -2.0 * g * b * alpha * alpha * mu_0 + 2.0 * g * alpha * mu_0 - 2.0 * alpha * lambda;
            let x8 = 2.0 * g * b * alpha * alpha * mu_0 + 4_f64 * g * alpha + 2.0 * alpha * lambda;
            let x9 = 2.0 * alpha + 2.0 * alpha * mu_0;
            let x10 = -2.0 * alpha;
            let x11 = 2.0 + mu_0;
            let x12 = -2.0 - mu_0;

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
    fn lim_coefficients_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;
        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);

        for i in 2..1000 {
            println!("-----");
            let alpha = PI / a * i as f64;

            let (
                (d_1_0, d_2_0, d_3_0, d_4_0, f_1_0, f_2_0, f_3_0, f_4_0),
                (d_1_1, d_2_1, d_3_1, d_4_1, f_1_1, f_2_1, f_3_1, f_4_1),
            ) = coefficients_1(b, alpha, mu_0, g, lambda);

            println!(
                "d_1_0: {d_1_0}, exp: {0}",
                (1.0 + mu_0) * 2.0 * g * b * mu_0 / (g + lambda) / (-2.0 - mu_0)
            );
            println!(
                "d_2_0: {d_2_0}, exp: {0}",
                (1.0 + mu_0) * 2.0 * b * mu_0 / (g + lambda) / (-2.0 - mu_0)
            );
            println!(
                "d_3_0: {d_3_0}, exp: {0}",
                -(1.0 + mu_0) * 2.0 * g * b * mu_0 / (g + lambda) / (-2.0 - mu_0)
            );
            println!(
                "d_4_0: {d_4_0}, exp: {0}",
                -(1.0 + mu_0) * 2.0 * b * mu_0 / (g + lambda) / (-2.0 - mu_0)
            );
            //
            println!(
                "f_1_0: {f_1_0}, exp: {0}",
                -(1.0 + mu_0) * 2.0 * g * b * mu_0 / (g + lambda) / (-2.0 - mu_0)
            );
            println!(
                "f_2_0: {f_2_0}, exp: {0}",
                -(1.0 + mu_0) * 2.0 * b * mu_0 / (g + lambda) / (-2.0 - mu_0)
            );
            println!(
                "f_3_0: {f_3_0}, exp: {0}",
                -(1.0 + mu_0) * 2.0 * g * b * mu_0 / (g + lambda) / (-2.0 - mu_0)
            );
            println!(
                "f_4_0: {f_4_0}, exp: {0}",
                -(1.0 + mu_0) * 2.0 * b * mu_0 / (g + lambda) / (-2.0 - mu_0)
            );
            //
            println!(
                "d_1_1: {d_1_1}, exp: {0}",
                4.0 * g * b * b * mu_0 * mu_0 / (g + lambda) / (-2.0 - mu_0)
            );
            println!(
                "d_2_1: {d_2_1}, exp: {0}",
                -8.0 * g * b * b * mu_0 * mu_0 / (g + lambda) / (-2.0 - mu_0)
            );
            println!(
                "d_3_1: {d_3_1}, exp: {0}",
                -4.0 * g * b * b * mu_0 * mu_0 / (g + lambda) / (-2.0 - mu_0)
            );
            println!(
                "d_4_1: {d_4_1}, exp: {0}",
                8.0 * g * b * b * mu_0 * mu_0 / (g + lambda) / (-2.0 - mu_0)
            );
            //
            println!("f_1_1: {f_1_1}, exp: {0}", 2.0);
            println!("f_2_1: {f_2_1}, exp: {0}", 4.0 / (2.0 + mu_0));
            println!(
                "f_3_1: {f_3_1}, exp: {0}",
                4.0 * g * b * b * mu_0 * mu_0 / (g + lambda) / (2.0 + mu_0)
            );
            println!(
                "f_4_1: {f_4_1}, exp: {0}",
                -4.0 * (1.0 + mu_0) / (2.0 + mu_0)
            );
        }
    }

    #[test]
    fn lim_psi_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;

        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);

        let x = 1.0;
        for i in 0..1111 {
            println!("----");
            let alpha = PI / a * i as f64;

            let ((psi_1_0, psi_2_0, psi_3_0, psi_4_0), (psi_1_1, psi_2_1, psi_3_1, psi_4_1)) =
                psi_1(x, b, alpha, mu_0, g, lambda);

            println!(
                "psi_1_0: {psi_1_0}, exp: {0}",
                -g * b * mu_0 / (g + lambda) / 2.0
            );
            println!(
                "psi_2_0: {psi_2_0}, exp: {0}",
                -b * mu_0 / (g + lambda) / 2.0
            );
            println!(
                "psi_3_0: {psi_3_0}, exp: {0}",
                -g * b * mu_0 / (g + lambda) / 2.0
            );
            println!(
                "psi_4_0: {psi_4_0}, exp: {0}",
                -b * mu_0 / (g + lambda) / 2.0
            );
            //
            println!(
                "psi_1_1: {psi_1_1}, exp: {0}",
                x * mu_0 / (1.0 + mu_0) / 2.0
            );
            println!("psi_2_1: {psi_2_1}, exp: {0}", x * mu_0 / (1.0 + mu_0));
            println!(
                "psi_3_1: {psi_3_1}, exp: {0}",
                -x * mu_0 / (1.0 + mu_0) / 2.0
            );
            println!("psi_4_1: {psi_4_1}, exp: {0}", x * mu_0 / (1.0 + mu_0));
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

        for x in 2..(b - 1_f64) as usize {
            let x = x as f64;
            for i in 2..20 {
                println!("----");
                let alpha = PI / a * i as f64;

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
            let alpha = PI / a * i as f64;

            let ((psi_1_0, psi_2_0, psi_3_0, psi_4_0), (psi_1_1, psi_2_1, psi_3_1, psi_4_1)) =
                psi_2(b, b, alpha, mu_0, g, lambda);
            let (
                (der_psi_1_0, der_psi_2_0, der_psi_3_0, der_psi_4_0),
                (der_psi_1_1, der_psi_2_1, der_psi_3_1, der_psi_4_1),
            ) = der_psi_2(b, b, alpha, mu_0, g, lambda);

            let line_1 = der_psi_1_0 + (-alpha) * psi_3_0;
            let line_2 = der_psi_2_0 + (-alpha) * psi_4_0;
            let line_3 = (2.0 * g + lambda) * der_psi_3_0 + (alpha * lambda) * psi_1_0;
            let line_4 = (2.0 * g + lambda) * der_psi_4_0 + (alpha * lambda) * psi_2_0;

            println!("0 line 1: {line_1}");
            println!("0 line 2: {line_2}");
            println!("0 line 3: {line_3}");
            println!("0 line 4: {line_4}");

            let line_1 = der_psi_1_1 + (-alpha) * psi_3_1;
            let line_2 = der_psi_2_1 + (-alpha) * psi_4_1;
            let line_3 = (2.0 * g + lambda) * der_psi_3_1 + (alpha * lambda) * psi_1_1;
            let line_4 = (2.0 * g + lambda) * der_psi_4_1 + (alpha * lambda) * psi_2_1;

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
            let alpha = PI / a * i as f64;

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

    #[test]
    fn a_1_test() {
        let a = 10_f64;
        let load_function = |x| x * x;
        let eps = 0.01;

        let x = 2.5;
        let a_1 = a_1(a, x, &load_function, eps);
        println!("{a_1}");
    }

    #[test]
    fn a_3_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;

        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);

        let x = 2.5;
        let ksi = 2.5;
        let a_3 = a_3(a, b, x, ksi, mu_0, g, lambda);
        println!("{a_3}");
    }

    #[test]
    fn f_m_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let mu_0 = mu_0(puasson_coef);
        let load_function = |x| x * x;
        let eps = 0.01;

        for m in 0..3 {
            let f_m = f_m(m, a, b, mu_0, &load_function, eps);
            println!("{f_m}");
        }
    }

    #[test]
    fn g_k_m_test() {
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
                let g_k_m = g_k_m(m, k, a, b, mu_0, g, lambda, eps);
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
        println!("{phi:?}");
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

        let h = 0.6;

        let unkn = unknown_function(h, a, b, mu_0, g, lambda, &load_function, eps);
        println!("{unkn}");
    }

    #[test]
    fn function_un_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;
        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);
        let load_function = |x| x * x;
        let eps = 0.1;

        let y = 5.0;
        let n = 1;
        let alpha = PI / a * n as f64;
        let function_un = function_un(a, b, y, alpha, mu_0, g, lambda, &load_function, eps);
        println!("{function_un}");
    }

    #[test]
    fn function_u_test() {
        let a = 10_f64;
        let b = 15_f64;
        let puasson_coef = 0.25;
        let young_modulus = 200_f64;
        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);
        let load_function = |x| x * x;
        let eps = 0.1;

        let x = 9.5;
        let y = 5.0;

        let function_u = function_u(a, b, x, y, mu_0, g, lambda, &load_function, eps);
        println!("{function_u}");
    }
}
