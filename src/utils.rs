use std::process::Command;

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
) -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
    let h_x = a / n_x as f64;
    let h_y = b / n_y as f64;

    let mut x = vec![0_f64; n_x as usize];
    let mut y = vec![0_f64; n_y as usize];
    let mut z = vec![vec![0_f64; n_x as usize]; n_y as usize];

    for j in 0..n_y {
        y[j as usize] = j as f64 * h_y;
        for i in 0..n_x {
            x[i as usize] = i as f64 * h_x;

            z[j as usize][i as usize] = f(x[i as usize], y[j as usize]);
            print!(
                "\r Calculating {}%...",
                f64::trunc(
                    (i + 1) as f64 * (j + 1) as f64 / (n_x as f64 * n_y as f64) * 100_f64 * 100_f64
                ) / 100_f64
            );
        }
    }
    (x, y, z)
}

pub fn surface_plot(x: Vec<f64>, y: Vec<f64>, z: Vec<Vec<f64>>) {
    Command::new("python3")
        .args([
            "plot/surface_plot.py",
            format!("-x={:?}", x).as_str(),
            format!("-y={:?}", y).as_str(),
            format!("-z={:?}", z).as_str(),
        ])
        .output()
        .unwrap();
}

#[test]
fn surface_plot_test() {
    surface_plot(
        vec![1_f64, 2_f64, 3_f64],
        vec![1_f64, 2_f64, 3_f64],
        vec![
            vec![5_f64, 5_f64, 5_f64],
            vec![5_f64, 5_f64, 5_f64],
            vec![5_f64, 5_f64, 5_f64],
        ],
    );
}
