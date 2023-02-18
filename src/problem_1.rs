use crate::{
    integration::definite_integral,
    utils::{function_calculation, g, lambda, mu_0, surface_static_plot},
    FunctionType, LoadFunction,
};
use clap::Parser;
use std::f64::consts::PI;

#[derive(Parser)]
#[clap(rename_all = "snake-case")]
pub struct Problem1 {
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

impl Problem1 {
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

        surface_static_plot(&x, &y, &z[0]);
    }
}

// returns c1, c2, c3, c4 coefficients
fn coefficients<F: Fn(f64) -> f64>(
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

fn function_un<F: Fn(f64) -> f64>(
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

fn function_vn<F: Fn(f64) -> f64>(
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

fn function_derivative_vn<F: Fn(f64) -> f64>(
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

    let res = c1 * e1 * (-y * alpha * alpha * mu_0 - alpha * mu_0)
        + c2 * e1 * (-y * alpha * alpha * mu_0 + 2_f64 * alpha)
        + c3 * e2 * (-y * alpha * alpha * mu_0 + alpha * mu_0)
        + c4 * e2 * (y * alpha * alpha * mu_0 + 2_f64 * alpha);

    res * coef
}

fn function_u<F: Fn(f64) -> f64>(
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
    let mut n = 10;
    let mut result = 0_f64;
    let mut prev_result;

    if x != a && x != 0_f64 {
        for i in 1..n {
            let alpha = PI * i as f64 / a;
            result += 2_f64
                * function_un(a, b, y, alpha, mu_0, g, lambda, load_function, eps)
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
                    * function_un(a, b, y, alpha, mu_0, g, lambda, load_function, eps)
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

fn function_derivative_u_x<F: Fn(f64) -> f64>(
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
    let mut n = 10;
    let mut result = 0_f64;
    let mut prev_result;

    if x != a && x != 0_f64 {
        for i in 1..n {
            let alpha = PI * i as f64 / a;
            result += 2_f64
                * alpha
                * function_un(a, b, y, alpha, mu_0, g, lambda, load_function, eps)
                * f64::cos(alpha * x)
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
                    * alpha
                    * function_un(a, b, y, alpha, mu_0, g, lambda, load_function, eps)
                    * f64::cos(alpha * x)
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

fn function_v<F: Fn(f64) -> f64>(
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
    let mut n = 10;
    let p0 = definite_integral(0_f64, a, 100, eps, load_function);
    let mut result = -p0 * y / a / (2_f64 * g + lambda);
    let mut prev_result;

    for i in 1..n {
        let alpha = PI * i as f64 / a;
        result += 2_f64
            * function_vn(a, b, y, alpha, mu_0, g, lambda, load_function, eps)
            * f64::cos(alpha * x)
            / a
            / (1_f64 + mu_0);
    }

    loop {
        n *= 2;
        prev_result = result;
        result = -p0 * y / a / (2_f64 * g + lambda);

        for i in 1..n {
            let alpha = PI * i as f64 / a;
            result += 2_f64
                * function_vn(a, b, y, alpha, mu_0, g, lambda, load_function, eps)
                * f64::cos(alpha * x)
                / a
                / (1_f64 + mu_0);
        }

        if f64::abs(result - prev_result) < eps {
            break;
        }
    }

    result
}

fn function_derivative_v_y<F: Fn(f64) -> f64>(
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
    let mut n = 10;
    let p0 = definite_integral(0_f64, a, 100, eps, load_function);
    let mut result = -p0 * y / a / (2_f64 * g + lambda);
    let mut prev_result;

    for i in 1..n {
        let alpha = PI * i as f64 / a;
        result += 2_f64
            * function_derivative_vn(a, b, y, alpha, mu_0, g, lambda, load_function, eps)
            * f64::cos(alpha * x)
            / a
            / (1_f64 + mu_0);
    }

    loop {
        n *= 2;
        prev_result = result;
        result = -p0 * y / a / (2_f64 * g + lambda);

        for i in 1..n {
            let alpha = PI * i as f64 / a;
            result += 2_f64
                * function_derivative_vn(a, b, y, alpha, mu_0, g, lambda, load_function, eps)
                * f64::cos(alpha * x)
                / a
                / (1_f64 + mu_0);
        }

        if f64::abs(result - prev_result) < eps {
            break;
        }
    }

    result
}

fn function_sigma_x<F: Fn(f64) -> f64>(
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

fn function_sigma_y<F: Fn(f64) -> f64>(
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
