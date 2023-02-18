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

pub fn function_calculation<F: Fn(f64, f64, f64) -> f64>(
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

    let mut x = vec![0_f64; n_x as usize];
    let mut y = vec![0_f64; n_y as usize];
    let mut t = vec![0_f64; n_t as usize];
    let mut z = vec![vec![vec![0_f64; n_x as usize]; n_y as usize]; n_t as usize];

    println!("Calculating ...");
    for k in 0..n_t {
        t[k as usize] = k as f64 * h_t;
        for j in 0..n_y {
            y[j as usize] = j as f64 * h_y;
            for i in 0..n_x {
                x[i as usize] = i as f64 * h_x;

                z[k as usize][j as usize][i as usize] =
                    f(x[i as usize], y[j as usize], t[k as usize]);
            }
        }
    }
    (x, y, z)
}

pub fn surface_static_plot(x: &[f64], y: &[f64], z: &[Vec<f64>]) {
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
