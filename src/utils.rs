use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use std::{env::current_dir, fs::File, io::Write, process::Command};

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
    initial_value: f64,
    f: &F,
    eps: f64,
    start: usize,
    mut n: usize,
) -> f64 {
    let mut result = initial_value;
    let mut prev_result;

    result += (start..n).into_par_iter().map(|i| f(i)).sum::<f64>();

    loop {
        n *= 2;
        prev_result = result;
        result = initial_value;

        result += (start..n).into_par_iter().map(|i| f(i)).sum::<f64>();

        if f64::abs(result - prev_result) < eps {
            break;
        }
    }
    result
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
    let x: Vec<_> = (0..n_x).into_par_iter().map(|i| i as f64 * h_x).collect();
    let y: Vec<_> = (0..n_y).into_par_iter().map(|i| i as f64 * h_y).collect();
    let t: Vec<_> = (0..n_t).into_par_iter().map(|i| i as f64 * h_t).collect();

    let z = t
        .into_par_iter()
        .map(|t| {
            y.clone()
                .into_par_iter()
                .map(|y| {
                    x.clone()
                        .into_par_iter()
                        .map(|x| {
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
    (x, y, z)
}

pub fn surface_static_plot(x: &[f64], y: &[f64], z: &[Vec<f64>]) {
    print!("Drawing ...");
    let path = current_dir().unwrap().join("result_tmp");
    let mut file = File::create(&path).unwrap();
    file.write_fmt(format_args!("{x:?}|{y:?}|{z:?}")).unwrap();

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

pub fn surface_dynamic_plot(x: &[f64], y: &[f64], z: &[Vec<Vec<f64>>]) {
    print!("Drawing ...");
    let path = current_dir().unwrap().join("result_tmp");
    let mut file = File::create(&path).unwrap();
    file.write_fmt(format_args!("{x:?}|{y:?}|{z:?}")).unwrap();

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
    use super::{surface_dynamic_plot, surface_static_plot};

    #[test]
    fn surface_static_plot_test() {
        surface_static_plot(
            &vec![1_f64, 2_f64, 3_f64],
            &vec![1_f64, 2_f64, 3_f64],
            &vec![
                vec![5_f64, 5_f64, 5_f64],
                vec![5_f64, 5_f64, 5_f64],
                vec![5_f64, 5_f64, 5_f64],
            ],
        );
    }

    #[test]
    fn surface_dynamic_plot_test() {
        surface_dynamic_plot(
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
        );
    }
}
