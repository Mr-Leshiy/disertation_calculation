// calculate Lamme's coefficient G value
pub fn g(puasson_coef: f64, young_modulus: f64) -> f64 {
    young_modulus / 2_f64 / (1_f64 + puasson_coef)
}

pub fn mu_0(puasson_coef: f64) -> f64 {
    1_f64 / (1_f64 - 2_f64 * puasson_coef)
}

pub fn lambda(puasson_coef: f64, young_modulus: f64) -> f64 {
    puasson_coef * young_modulus / (1_f64 + puasson_coef) / (1_f64 - 2_f64 * puasson_coef)
}

pub fn function_calculation<F: Fn(f64, f64) -> f64>(
    a: f64,
    b: f64,
    n_x: u32,
    n_y: u32,
    f: F,
) -> Vec<Vec<f64>> {
    let h_x = a / n_x as f64;
    let h_y = b / n_y as f64;

    let mut res = vec![vec![0_f64; n_x as usize]; n_y as usize];

    for j in 0..n_y {
        for i in 0..n_x {
            let x = i as f64 * h_x;
            let y = j as f64 * h_y;

            res[j as usize][i as usize] = f(x, y);
            print!(
                "\r Calculating {}%...",
                (i + 1) as f64 * (j + 1) as f64 / (n_x as f64 * n_y as f64) * 100_f64
            );
        }
    }
    res
}
