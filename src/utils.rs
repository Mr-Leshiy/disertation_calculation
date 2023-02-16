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
