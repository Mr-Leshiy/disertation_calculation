use rayon::prelude::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use std::{env::current_dir, fs::File, io::Write, path::PathBuf, process::Command};

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

pub fn sum_calc<F: Fn(usize) -> f64 + Send + Sync>(
    f: &F,
    eps: f64,
    start: usize,
    mut n: usize,
) -> f64 {
    let mut result = 0.0;
    let mut prev_result;

    result += (start..n).into_par_iter().map(|i| f(i)).sum::<f64>();

    loop {
        n *= 2;
        prev_result = result;
        result = 0.0;

        result += (start..n).into_par_iter().map(|i| f(i)).sum::<f64>();

        if f64::abs(result - prev_result) < eps {
            break;
        }
    }
    result
}

pub fn sum_calc_finit<F: Fn(usize) -> f64 + Send + Sync>(f: &F, start: usize, n: usize) -> f64 {
    (start..n).into_par_iter().map(|i| f(i)).sum::<f64>()
}

pub fn function_calculation<F: Fn(f64, f64, f64) -> f64 + Send + Sync>(
    a: f64,
    b: f64,
    n_x: u32,
    n_y: u32,
    t_data: Option<(f64, u32)>,
    f: F,
) -> (Vec<f64>, Vec<f64>, Vec<Vec<Vec<f64>>>) {
    let (t, n_t) = t_data.unwrap_or((1_f64, 1_u32));

    let h_x = a / n_x as f64;
    let h_y = b / n_y as f64;
    let h_t = t / n_t as f64;

    println!("Calculating ...\n");
    let x: Vec<_> = (0..n_x + 1)
        .into_par_iter()
        .map(|i| i as f64 * h_x)
        .collect();
    let y: Vec<_> = (0..n_y + 1)
        .into_par_iter()
        .map(|i| i as f64 * h_y)
        .collect();
    let t: Vec<_> = (0..n_t + 1)
        .into_par_iter()
        .map(|i| i as f64 * h_t)
        .collect();

    let z = t
        .into_par_iter()
        .map(|t| {
            y.clone()
                .into_par_iter()
                .map(|y| {
                    x.clone()
                        .into_par_iter()
                        .enumerate()
                        .map(|(_, x)| {
                            println!("{x}|{y}|{t}| - start");
                            let res = f(x, y, t);
                            println!("{x}|{y}|{t}| - done");
                            res
                        })
                        .collect()
                })
                .collect()
        })
        .collect();
    println!("Finished calculation");
    (x, y, z)
}

pub fn function_calculation2<F: Fn(f64) -> f64 + Send + Sync>(
    omega_1: f64,
    omega_2: f64,
    n_omega: u32,
    f: F,
) -> (Vec<f64>, Vec<f64>) {
    let h_omega = (omega_2 - omega_1) / n_omega as f64;

    println!("Calculating ...\n");
    let omega: Vec<_> = (0..n_omega + 1)
        .into_par_iter()
        .map(|i| i as f64 * h_omega + omega_1)
        .collect();

    let z = omega
        .clone()
        .into_iter()
        .map(|omega| {
            println!("{omega} - start");
            let res = f(omega);
            println!("{omega} - done");
            res
        })
        .collect();
    println!("Finished calculation");
    (omega, z)
}

pub fn save_freq(omega: &[f64], z: &[f64], file_name: &str) -> PathBuf {
    let path = current_dir().unwrap().join(file_name);
    let mut file = File::create(&path).unwrap();
    file.write_fmt(format_args!(
        "{}|{}",
        serde_json::to_string(&omega).unwrap(),
        serde_json::to_string(&z).unwrap()
    ))
    .unwrap();
    path
}

pub fn save_static(x: &[f64], y: &[f64], z: &[Vec<f64>], file_name: &str) -> PathBuf {
    let path = current_dir().unwrap().join(file_name);
    let mut file = File::create(&path).unwrap();
    file.write_fmt(format_args!(
        "{}|{}|{}",
        serde_json::to_string(&x).unwrap(),
        serde_json::to_string(&y).unwrap(),
        serde_json::to_string(&z).unwrap()
    ))
    .unwrap();
    path
}

pub fn save_dynamic(x: &[f64], y: &[f64], z: &[Vec<Vec<f64>>], file_name: &str) -> PathBuf {
    let path = current_dir().unwrap().join(file_name);
    let mut file = File::create(&path).unwrap();
    file.write_fmt(format_args!(
        "{}|{}|{}",
        serde_json::to_string(&x).unwrap(),
        serde_json::to_string(&y).unwrap(),
        serde_json::to_string(&z).unwrap()
    ))
    .unwrap();
    path
}

pub fn function_plot(path: &PathBuf) {
    print!("Drawing ...");

    let _res = Command::new("python3")
        .args([
            "plot/function_plot.py",
            format!("-p={}", path.to_str().unwrap()).as_str(),
        ])
        .output()
        .unwrap();
    println!(
        "stout: {}, stderr: {}",
        String::from_utf8(_res.stdout).unwrap(),
        String::from_utf8(_res.stderr).unwrap()
    );
}

pub fn surface_static_plot(path: &PathBuf) {
    print!("Drawing ...");

    let _res = Command::new("python3")
        .args([
            "plot/surface_static_plot.py",
            format!("-p={}", path.to_str().unwrap()).as_str(),
        ])
        .output()
        .unwrap();
    println!(
        "stout: {}, stderr: {}",
        String::from_utf8(_res.stdout).unwrap(),
        String::from_utf8(_res.stderr).unwrap()
    );
}

pub fn surface_dynamic_plot(path: &PathBuf) {
    print!("Drawing ...");

    let _res = Command::new("python3")
        .args([
            "plot/surface_dynamic_plot.py",
            format!("-p={}", path.to_str().unwrap()).as_str(),
        ])
        .output()
        .unwrap();
    println!(
        "stout: {}, stderr: {}",
        String::from_utf8(_res.stdout).unwrap(),
        String::from_utf8(_res.stderr).unwrap()
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn function_plot_test() {
        let path = save_freq(
            &vec![1_f64, 2_f64, 3_f64],
            &vec![1_f64, 2_f64, 3_f64],
            "tmp",
        );
        function_plot(&path);
    }

    #[test]
    fn surface_static_plot_test() {
        let path = save_static(
            &vec![1_f64, 2_f64, 3_f64],
            &vec![1_f64, 2_f64, 3_f64],
            &vec![
                vec![5_f64, 5_f64, 5_f64],
                vec![5_f64, 5_f64, 5_f64],
                vec![5_f64, 5_f64, 5_f64],
            ],
            "tmp",
        );
        surface_static_plot(&path);
    }

    #[test]
    fn surface_dynamic_plot_test() {
        let path = save_dynamic(
            &vec![1_f64, 2_f64, 3_f64],
            &vec![1_f64, 2_f64, 3_f64],
            &vec![
                vec![
                    vec![5_f64, 5_f64, 5_f64],
                    vec![5_f64, 5_f64, 5_f64],
                    vec![5_f64, 5_f64, 5_f64],
                ],
                vec![
                    vec![6_f64, 6_f64, 6_f64],
                    vec![6_f64, 6_f64, 6_f64],
                    vec![6_f64, 6_f64, 6_f64],
                ],
            ],
            "tmp",
        );

        surface_dynamic_plot(&path);
    }
}
