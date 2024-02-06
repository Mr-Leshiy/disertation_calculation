#![allow(clippy::too_many_arguments)]
#![allow(dead_code)]

use utils::{function_calculation, g, lambda, mu_0, FunctionData, FunctionFileInfo};

mod dynamic_1;
mod integration;
mod matrices;
mod polynomials;
mod static_1;
mod static_2;
mod utils;

trait LoadFuncT: Fn(f64) -> f64 + Send + Sync {}
impl<F: Fn(f64) -> f64 + Send + Sync> LoadFuncT for F {}

struct FunctionCalculation<Func: LoadFuncT> {
    a1: f64,
    a2: f64,

    b1: f64,
    b2: f64,

    puasson_coef: f64,
    young_modulus: f64,

    g: f64,
    lambda: f64,
    mu_0: f64,

    load_function: Func,

    eps: f64,
}

impl<Func: LoadFuncT> FunctionCalculation<Func> {
    fn new(
        a1: f64,
        a2: f64,
        b1: f64,
        b2: f64,
        puasson_coef: f64,
        young_modulus: f64,
        load_function: Func,
        eps: f64,
    ) -> Self {
        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);

        Self {
            a1,
            a2,
            b1,
            b2,
            puasson_coef,
            young_modulus,
            g,
            lambda,
            mu_0,
            load_function,
            eps,
        }
    }

    fn static_1_function_u_fixed_x(&self, x: &[f64]) -> FunctionFileInfo {
        let mut data = Vec::new();
        for x in x {
            let (res_x, res_y) = function_calculation(self.b1, self.b2, 100, |y| {
                static_1::function_u(
                    self.a2,
                    self.b2,
                    *x,
                    y,
                    self.mu_0,
                    self.g,
                    self.lambda,
                    &self.load_function,
                    self.eps,
                )
            });
            data.push(FunctionData {
                x: res_x,
                y: res_y,
                label: format!("a = {}; b = {}, x = {}", self.a2, self.b2, x),
            });
        }

        FunctionFileInfo {
            data,
            ox_name: "y".to_string(),
            oy_name: "u(x,y)".to_string(),
            file_name: format!("static_1_function_u_x"),
        }
    }

    fn static_1_function_u_fixed_y(&self, y: &[f64]) -> FunctionFileInfo {
        let mut data = Vec::new();
        for y in y {
            let (res_x, res_y) = function_calculation(self.a1, self.a2, 100, |x| {
                static_1::function_u(
                    self.a2,
                    self.b2,
                    x,
                    *y,
                    self.mu_0,
                    self.g,
                    self.lambda,
                    &self.load_function,
                    self.eps,
                )
            });
            data.push(FunctionData {
                x: res_x,
                y: res_y,
                label: format!("a = {}; b = {}, y = {}", self.a2, self.b2, y),
            });
        }

        FunctionFileInfo {
            data,
            ox_name: "x".to_string(),
            oy_name: "u(x,y)".to_string(),
            file_name: format!("static_1_function_u_y"),
        }
    }

    fn static_1_function_v_fixed_x(&self, x: &[f64]) -> FunctionFileInfo {
        let mut data = Vec::new();
        for x in x {
            let (res_x, res_y) = function_calculation(self.b1, self.b2, 100, |y| {
                static_1::function_v(
                    self.a2,
                    self.b2,
                    *x,
                    y,
                    self.mu_0,
                    self.g,
                    self.lambda,
                    &self.load_function,
                    self.eps,
                )
            });
            data.push(FunctionData {
                x: res_x,
                y: res_y,
                label: format!("a = {}; b = {}, x = {}", self.a2, self.b2, x),
            });
        }

        FunctionFileInfo {
            data,
            ox_name: "y".to_string(),
            oy_name: "v(x,y)".to_string(),
            file_name: format!("static_1_function_u_x"),
        }
    }

    fn static_1_function_v_fixed_y(&self, y: &[f64]) -> FunctionFileInfo {
        let mut data = Vec::new();
        for y in y {
            let (res_x, res_y) = function_calculation(self.a1, self.a2, 100, |x| {
                static_1::function_v(
                    self.a2,
                    self.b2,
                    x,
                    *y,
                    self.mu_0,
                    self.g,
                    self.lambda,
                    &self.load_function,
                    self.eps,
                )
            });
            data.push(FunctionData {
                x: res_x,
                y: res_y,
                label: format!("a = {}; b = {}, y = {}", self.a2, self.b2, y),
            });
        }

        FunctionFileInfo {
            data,
            ox_name: "x".to_string(),
            oy_name: "v(x,y)".to_string(),
            file_name: format!("static_1_function_u_y"),
        }
    }

    fn static_1_function_tau_xy_fixed_x(&self, x: &[f64]) -> FunctionFileInfo {
        let mut data = Vec::new();
        for x in x {
            let (res_x, res_y) = function_calculation(self.b1, self.b2, 100, |y| {
                static_1::function_tau_xy(
                    self.a2,
                    self.b2,
                    *x,
                    y,
                    self.mu_0,
                    self.g,
                    self.lambda,
                    &self.load_function,
                    self.eps,
                )
            });
            data.push(FunctionData {
                x: res_x,
                y: res_y,
                label: format!("a = {}; b = {}, x = {}", self.a2, self.b2, x),
            });
        }

        FunctionFileInfo {
            data,
            ox_name: "y".to_string(),
            oy_name: "tau_xy(x,y)".to_string(),
            file_name: format!("static_1_function_tau_x"),
        }
    }

    fn static_1_function_tau_xy_fixed_y(&self, y: &[f64]) -> FunctionFileInfo {
        let mut data = Vec::new();
        for y in y {
            let (res_x, res_y) = function_calculation(self.a1, self.a2, 100, |x| {
                static_1::function_tau_xy(
                    self.a2,
                    self.b2,
                    x,
                    *y,
                    self.mu_0,
                    self.g,
                    self.lambda,
                    &self.load_function,
                    self.eps,
                )
            });
            data.push(FunctionData {
                x: res_x,
                y: res_y,
                label: format!("a = {}; b = {}, y = {}", self.a2, self.b2, y),
            });
        }

        FunctionFileInfo {
            data,
            ox_name: "x".to_string(),
            oy_name: "tau_xy(x,y)".to_string(),
            file_name: format!("static_1_function_tau_y"),
        }
    }

    fn static_1_function_sigma_x_fixed_x(&self, x: &[f64]) -> FunctionFileInfo {
        let mut data = Vec::new();
        for x in x {
            let (res_x, res_y) = function_calculation(self.b1, self.b2, 100, |y| {
                static_1::function_sigma_x(
                    self.a2,
                    self.b2,
                    *x,
                    y,
                    self.mu_0,
                    self.g,
                    self.lambda,
                    &self.load_function,
                    self.eps,
                )
            });
            data.push(FunctionData {
                x: res_x,
                y: res_y,
                label: format!("a = {}; b = {}, x = {}", self.a2, self.b2, x),
            });
        }

        FunctionFileInfo {
            data,
            ox_name: "y".to_string(),
            oy_name: "sigma_x(x,y)".to_string(),
            file_name: format!("static_1_function_sigma_x_x"),
        }
    }

    fn static_1_function_sigma_x_fixed_y(&self, y: &[f64]) -> FunctionFileInfo {
        let mut data = Vec::new();
        for y in y {
            let (res_x, res_y) = function_calculation(self.a1, self.a2, 100, |x| {
                static_1::function_sigma_x(
                    self.a2,
                    self.b2,
                    x,
                    *y,
                    self.mu_0,
                    self.g,
                    self.lambda,
                    &self.load_function,
                    self.eps,
                )
            });
            data.push(FunctionData {
                x: res_x,
                y: res_y,
                label: format!("a = {}; b = {}, y = {}", self.a2, self.b2, y),
            });
        }

        FunctionFileInfo {
            data,
            ox_name: "x".to_string(),
            oy_name: "sigma_x(x,y)".to_string(),
            file_name: format!("static_1_function_sigma_x_y"),
        }
    }


    fn static_1_function_sigma_y_fixed_x(&self, x: &[f64]) -> FunctionFileInfo {
        let mut data = Vec::new();
        for x in x {
            let (res_x, res_y) = function_calculation(self.b1, self.b2, 100, |y| {
                static_1::function_sigma_y(
                    self.a2,
                    self.b2,
                    *x,
                    y,
                    self.mu_0,
                    self.g,
                    self.lambda,
                    &self.load_function,
                    self.eps,
                )
            });
            data.push(FunctionData {
                x: res_x,
                y: res_y,
                label: format!("a = {}; b = {}, x = {}", self.a2, self.b2, x),
            });
        }

        FunctionFileInfo {
            data,
            ox_name: "y".to_string(),
            oy_name: "sigma_y(x,y)".to_string(),
            file_name: format!("static_1_function_sigma_y_x"),
        }
    }

    fn static_1_function_sigma_y_fixed_y(&self, y: &[f64]) -> FunctionFileInfo {
        let mut data = Vec::new();
        for y in y {
            let (res_x, res_y) = function_calculation(self.a1, self.a2, 100, |x| {
                static_1::function_sigma_y(
                    self.a2,
                    self.b2,
                    x,
                    *y,
                    self.mu_0,
                    self.g,
                    self.lambda,
                    &self.load_function,
                    self.eps,
                )
            });
            data.push(FunctionData {
                x: res_x,
                y: res_y,
                label: format!("a = {}; b = {}, y = {}", self.a2, self.b2, y),
            });
        }

        FunctionFileInfo {
            data,
            ox_name: "x".to_string(),
            oy_name: "sigma_y(x,y)".to_string(),
            file_name: format!("static_1_function_sigma_y_y"),
        }
    }
}

#[cfg(test)]
mod run {
    use super::*;
    use crate::utils::{function_plot, save_function};

    #[test]
    fn run() {
        let b = 15.0;
        let a = b / 2.0;

        // steel
        let puasson_coef = 0.25;
        let young_modulus = 200.0;

        let load_function = |x| (x - 5.5) * (x - 5.5);
        let eps = 0.1;

        let func_calc = FunctionCalculation::new(
            0.0,
            a,
            0.0,
            b,
            puasson_coef,
            young_modulus,
            load_function,
            eps,
        );

        let res = func_calc.static_1_function_sigma_y_fixed_y(&[b]);

        let file = save_function(res);
        function_plot(&file)
    }
}
