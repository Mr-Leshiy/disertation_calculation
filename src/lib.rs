#![allow(clippy::too_many_arguments)]
#![allow(dead_code)]

mod dynamic_1;
mod integration;
mod matrices;
mod polynomials;
mod static_1;
mod static_2;
mod utils;

#[cfg(test)]
mod run {
    use super::*;
    use crate::utils::{
        function_calculation, function_plot, g, lambda, mu_0, save_function, FunctionData,
    };

    #[test]
    fn run() {
        let a = 10.0;
        let b = 15.0;
        // steel
        let puasson_coef = 0.25;
        let young_modulus = 200.0;

        let g = g(puasson_coef, young_modulus);
        let lambda = lambda(puasson_coef, young_modulus);
        let mu_0 = mu_0(puasson_coef);

        let load_function = |x| (x - 5.5) * (x - 5.5);
        let eps = 0.1;

        let n = 5;

        let y = 1.0;

        let (res_x, res_y) = function_calculation(-a, a, 100, |x| {
            static_2::function_u(a, b, x, y, mu_0, g, lambda, &load_function, eps)
        });
        let case1 = FunctionData {
            x: &res_x,
            y: &res_y,
            label: "",
        };

        let file_name = "test";
        let file = save_function(vec![case1], "x", "u(x,y)", file_name);
        function_plot(&file)
    }
}
