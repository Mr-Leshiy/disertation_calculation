use std::f64::consts::PI;

use clap::Parser;

use crate::{
    integration::definite_integral,
    utils::{function_calculation, g, lambda, mu_0, surface_dynamic_plot},
    FunctionType, LoadFunction,
};

#[derive(Parser)]
#[clap(rename_all = "snake-case")]
pub struct Problem3 {
    #[clap(long)]
    a: f64,
    #[clap(long)]
    b: f64,
    #[clap(long)]
    t: f64,
    #[clap(long)]
    n_x: u32,
    #[clap(long)]
    n_y: u32,
    #[clap(long)]
    n_t: u32,
    #[clap(long)]
    omega: f64,
    #[clap(long)]
    c1: f64,
    #[clap(long)]
    c2: f64,
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

impl Problem3 {
    pub fn exec(self) {
        assert!(self.a < self.b);

        let mu_0 = mu_0(self.puasson_coef);
        let g = g(self.puasson_coef, self.young_modulus);
        let lambda = lambda(self.puasson_coef, self.young_modulus);

        let (x, y, z) = match self.function_type {
            FunctionType::U => function_calculation(
                self.a,
                self.b,
                self.n_x,
                self.n_y,
                Some((self.t, self.n_t)),
                |x, y, t| {
                    function_u(
                        self.a,
                        self.b,
                        x,
                        y,
                        t,
                        self.omega,
                        self.c1,
                        self.c2,
                        mu_0,
                        g,
                        lambda,
                        &|x| self.load_function.call(x),
                        self.eps,
                    )
                },
            ),
            FunctionType::V => function_calculation(
                self.a,
                self.b,
                self.n_x,
                self.n_y,
                Some((self.t, self.n_t)),
                |x, y, t| {
                    function_v(
                        self.a,
                        self.b,
                        x,
                        y,
                        t,
                        self.omega,
                        self.c1,
                        self.c2,
                        mu_0,
                        g,
                        lambda,
                        &|x| self.load_function.call(x),
                        self.eps,
                    )
                },
            ),
            FunctionType::SigmaX => function_calculation(
                self.a,
                self.b,
                self.n_x,
                self.n_y,
                Some((self.t, self.n_t)),
                |x, y, t| {
                    function_sigma_x(
                        self.a,
                        self.b,
                        x,
                        y,
                        t,
                        self.omega,
                        self.c1,
                        self.c2,
                        mu_0,
                        g,
                        lambda,
                        &|x| self.load_function.call(x),
                        self.eps,
                    )
                },
            ),
            FunctionType::SigmaY => function_calculation(
                self.a,
                self.b,
                self.n_x,
                self.n_y,
                Some((self.t, self.n_t)),
                |x, y, t| {
                    function_sigma_y(
                        self.a,
                        self.b,
                        x,
                        y,
                        t,
                        self.omega,
                        self.c1,
                        self.c2,
                        mu_0,
                        g,
                        lambda,
                        &|x| self.load_function.call(x),
                        self.eps,
                    )
                },
            ),
        };

        surface_dynamic_plot(&x, &y, &z);
    }
}

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
        / (2_f64 * a1 * (a1 * a1 - a3 * a3))
        + c_2 * a1 * alpha * mu_0 * (e1 + e2) / (2_f64 * a1 * (a1 * a1 - a3 * a3))
        + c_3 * (a3 * a3 + a3 * a3 * mu_0 - alpha * alpha + omega * omega / c2 / c2) * (e3 - e4)
            / (2_f64 * a3 * (a3 * a3 - a1 * a1))
        + c_4 * a3 * alpha * mu_0 * (e3 + e4) / (2_f64 * a3 * (a3 * a3 - a1 * a1));

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

    let res = c_1 * (-a1 * alpha * mu_0) * (e1 + e2) / (2_f64 * a1 * (a1 * a1 - a3 * a3))
        + c_2
            * (a1 * a1 - alpha * alpha - alpha * alpha * mu_0 + omega * omega / c1 / c1)
            * (e1 - e2)
            / (2_f64 * a1 * (a1 * a1 - a3 * a3))
        + c_3 * (-a3 * alpha * mu_0) * (e3 + e4) / (2_f64 * a3 * (a3 * a3 - a1 * a1))
        + c_4
            * (a3 * a3 - alpha * alpha - alpha * alpha * mu_0 + omega * omega / c1 / c1)
            * (e3 - e4)
            / (2_f64 * a3 * (a3 * a3 - a1 * a1));

    res
}

fn function_derivative_vn<F: Fn(f64) -> f64>(
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

    let res = c_1 * (-a1 * alpha * mu_0) * a1 * (e1 - e2) / (2_f64 * a1 * (a1 * a1 - a3 * a3))
        + c_2
            * (a1 * a1 - alpha * alpha - alpha * alpha * mu_0 + omega * omega / c1 / c1)
            * a1
            * (e1 + e2)
            / (2_f64 * a1 * (a1 * a1 - a3 * a3))
        + c_3 * (-a3 * alpha * mu_0) * a3 * (e3 - e4) / (2_f64 * a3 * (a3 * a3 - a1 * a1))
        + c_4
            * (a3 * a3 - alpha * alpha - alpha * alpha * mu_0 + omega * omega / c1 / c1)
            * a3
            * (e3 + e4)
            / (2_f64 * a3 * (a3 * a3 - a1 * a1));

    res
}

fn function_u<F: Fn(f64) -> f64>(
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
    let mut n = 10;
    let mut result = 0_f64;
    let mut prev_result;

    // actually we should have a e^(i*omega*t), but for simplicity we are taking only real part
    let cos_t = f64::cos(omega * t);

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
                * cos_t
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
                    * cos_t
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
    let mut n = 10;
    let mut result = 0_f64;
    let mut prev_result;

    // actually we should have a e^(i*omega*t), but for simplicity we are taking only real part
    let cos_t = f64::cos(omega * t);

    if x != a && x != 0_f64 {
        for i in 1..n {
            let alpha = PI * i as f64 / a;
            result += 2_f64
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
                * cos_t
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
                    * cos_t
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
    let mut n = 10;
    let p0 = definite_integral(0_f64, a, 100, eps, load_function);
    let v0 = -p0 * f64::sin(y * f64::sqrt(omega * omega / (c2 * c2 * (1_f64 + mu_0))))
        / ((2_f64 * g + lambda)
            * f64::sqrt(omega * omega / (c2 * c2 * (1_f64 + mu_0)))
            * f64::sin(b * f64::sqrt(omega * omega / (c2 * c2 * (1_f64 + mu_0)))));
    let mut result = v0 / a;
    let mut prev_result;

    // actually we should have a e^(i*omega*t), but for simplicity we are taking only real part
    let cos_t = f64::cos(omega * t);

    for i in 1..n {
        let alpha = PI * i as f64 / a;
        result += 2_f64
            * function_vn(
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
            * cos_t
            / a
            / (1_f64 + mu_0);
    }

    loop {
        n *= 2;
        prev_result = result;
        result = v0 / a;

        for i in 1..n {
            let alpha = PI * i as f64 / a;
            result += 2_f64
                * function_vn(
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
                * cos_t
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
    let mut n = 10;
    let p0 = definite_integral(0_f64, a, 100, eps, load_function);
    let v0 = -p0 * f64::sin(y * f64::sqrt(omega * omega / (c2 * c2 * (1_f64 + mu_0))))
        / ((2_f64 * g + lambda)
            * f64::sqrt(omega * omega / (c2 * c2 * (1_f64 + mu_0)))
            * f64::sin(b * f64::sqrt(omega * omega / (c2 * c2 * (1_f64 + mu_0)))));
    let mut result = v0 / a;
    let mut prev_result;

    // actually we should have a e^(i*omega*t), but for simplicity we are taking only real part
    let cos_t = f64::cos(omega * t);

    for i in 1..n {
        let alpha = PI * i as f64 / a;
        result += 2_f64
            * function_derivative_vn(
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
            * cos_t
            / a
            / (1_f64 + mu_0);
    }

    loop {
        n *= 2;
        prev_result = result;
        result = v0 / a;

        for i in 1..n {
            let alpha = PI * i as f64 / a;
            result += 2_f64
                * function_derivative_vn(
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
                * cos_t
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

fn function_sigma_y<F: Fn(f64) -> f64>(
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

    #[test]
    fn a_coefficients_test() {
        let a = 10_f64;
        let omega = 0.75;
        let c1 = 10_f64;
        let c2 = 10_f64;
        let puasson_coef = 0.25;
        let eps = 0.001;

        let mu_0 = mu_0(puasson_coef);

        for i in 1..10 {
            let alpha = PI * i as f64 / a;

            let (a1, a3) = a_coefficients(omega, c1, c2, alpha, mu_0);

            let func = |s: f64| {
                (s * s - alpha * alpha - alpha * alpha * mu_0 + omega * omega / c1 / c1)
                    * (s * s + s * s * mu_0 - alpha * alpha + omega * omega / c2 / c2)
                    + s * s * alpha * alpha * mu_0 * mu_0
            };

            assert!(
                f64::abs(func(a1) - 0_f64) < eps,
                "result: {}, exp: {}",
                func(a1),
                0_f64
            );
            assert!(
                f64::abs(func(-a1) - 0_f64) < eps,
                "result: {}, exp: {}",
                func(-a1),
                0_f64
            );
            assert!(
                f64::abs(func(a3) - 0_f64) < eps,
                "result: {}, exp: {}",
                func(a3),
                0_f64
            );
            assert!(
                f64::abs(func(-a3) - 0_f64) < eps,
                "result: {}, exp: {}",
                func(-a3),
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

            let (a1, a3) = a_coefficients(omega, c1, c2, alpha, mu_0);

            let pn = definite_integral(0_f64, a, 100, eps, &|x| {
                load_function(x) * f64::cos(alpha * x)
            });

            let e1 = f64::exp(a1 * b) + f64::exp(-a1 * b);
            let e2 = f64::exp(a1 * b) - f64::exp(-a1 * b);
            let e3 = f64::exp(a3 * b) + f64::exp(-a3 * b);
            let e4 = f64::exp(a3 * b) - f64::exp(-a3 * b);

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
                f64::abs(eq1 - 0_f64) < eps,
                "result: {}, exp: {}",
                eq1,
                0_f64
            );
            assert!(
                f64::abs(eq2 + pn * (1_f64 + mu_0)) < eps,
                "result: {}, exp: {}",
                eq2,
                -pn * (1_f64 + mu_0)
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
