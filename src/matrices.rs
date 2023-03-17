use rayon::prelude::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

#[derive(Debug, Clone, PartialEq)]
pub struct Matrix {
    elements: Vec<Vec<f64>>,
}

impl Matrix {
    pub fn new(elements: Vec<Vec<f64>>) -> Self {
        let n = elements.len();
        assert!(n > 0);
        let m = elements[0].len();
        assert!(m > 0);
        elements.par_iter().for_each(|row| assert_eq!(row.len(), m));

        Self { elements }
    }

    pub fn n(&self) -> usize {
        self.elements.len()
    }

    pub fn m(&self) -> usize {
        self.elements[0].len()
    }

    pub fn get_element(&self, row: usize, column: usize) -> f64 {
        self.elements[row][column]
    }

    pub fn update_element(&mut self, row: usize, column: usize, val: f64) {
        self.elements[row][column] = val;
    }

    pub fn raw_swap(&mut self, row1: usize, row2: usize) {
        if row1 != row2 {
            (0..self.m()).into_iter().for_each(|column| {
                let tmp = self.elements[row1][column];
                self.elements[row1][column] = self.elements[row2][column];
                self.elements[row2][column] = tmp;
            });
        }
    }

    pub fn column_swap(&mut self, column1: usize, column2: usize) {
        if column1 != column2 {
            (0..self.n()).into_iter().for_each(|row| {
                let tmp = self.elements[row][column1];
                self.elements[row][column1] = self.elements[row][column2];
                self.elements[row][column2] = tmp;
            });
        }
    }
}

pub fn add(a: &Matrix, b: &Matrix) -> Matrix {
    assert_eq!(a.n(), b.n());
    assert_eq!(a.m(), b.m());

    let res = (0..a.n())
        .into_par_iter()
        .map(|row| {
            (0..a.m())
                .into_par_iter()
                .map(|column| a.get_element(row, column) + b.get_element(row, column))
                .collect()
        })
        .collect();

    Matrix::new(res)
}

pub fn sub(a: &Matrix, b: &Matrix) -> Matrix {
    assert_eq!(a.n(), b.n());
    assert_eq!(a.m(), b.m());

    let res = (0..a.n())
        .into_par_iter()
        .map(|row| {
            (0..b.m())
                .into_par_iter()
                .map(|column| a.get_element(row, column) - b.get_element(row, column))
                .collect()
        })
        .collect();

    Matrix::new(res)
}

pub fn mul(a: &Matrix, b: &Matrix) -> Matrix {
    assert_eq!(a.m(), b.n());

    let res: Vec<Vec<_>> = (0..a.n())
        .into_par_iter()
        .map(|row_1| {
            (0..b.m())
                .into_par_iter()
                .map(|column_2| {
                    (0..a.m())
                        .zip(0..b.n())
                        .map(|(comunm_1, row_2)| {
                            a.get_element(row_1, comunm_1) * b.get_element(row_2, column_2)
                        })
                        .sum::<f64>()
                })
                .collect()
        })
        .collect();

    Matrix::new(res)
}

pub fn gaussian_elimination(mut a: Matrix, mut left: Matrix, eps: f64) -> (Matrix, Matrix) {
    assert_eq!(left.n(), 1);
    assert_eq!(left.m(), a.n());

    let mut h = 0;
    let mut k = 0;
    while h < a.n() && k < a.m() {
        let mut max = f64::abs(a.get_element(h, k));
        let mut i_max = h;
        (h..a.n()).into_iter().for_each(|i| {
            if a.get_element(i, k) > max {
                max = a.get_element(i, k);
                i_max = i;
            }
        });
        if f64::abs(a.get_element(i_max, k)) < eps {
            k += 1;
        } else {
            a.raw_swap(h, i_max);
            left.column_swap(h, i_max);
            (h + 1..a.n()).into_iter().for_each(|i| {
                let f = a.get_element(i, k) / a.get_element(h, k);
                a.update_element(i, k, 0_f64);

                let left_new_val = left.get_element(0, i) - left.get_element(0, h) * f;
                left.update_element(0, i, left_new_val);

                (k + 1..a.m()).into_iter().for_each(|j| {
                    let new_val = a.get_element(i, j) - a.get_element(h, j) * f;
                    a.update_element(i, j, new_val);
                })
            });
            h += 1;
            k += 1;
        }
    }
    (a, left)
}

pub fn system_solve(a: Matrix, left: Matrix, eps: f64) -> Matrix {
    assert_eq!(a.n(), a.m());
    let (a, left) = gaussian_elimination(a, left, eps);
    let mut res = Matrix::new(vec![vec![0_f64; a.n()]]);

    (0..a.n()).rev().for_each(|i| {
        let sum = (i + 1..a.m())
            .rev()
            .map(|j| res.get_element(0, j) * a.get_element(i, j))
            .sum::<f64>();
        let left_val = left.get_element(0, i) - sum;
        let val = left_val / a.get_element(i, i);
        res.update_element(0, i, val);
    });
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matrix_raw_swap_test() {
        let mut a = Matrix::new(vec![
            vec![1_f64, 2_f64, 3_f64],
            vec![4_f64, 5_f64, 6_f64],
            vec![7_f64, 8_f64, 9_f64],
        ]);
        let expected = Matrix::new(vec![
            vec![7_f64, 8_f64, 9_f64],
            vec![4_f64, 5_f64, 6_f64],
            vec![1_f64, 2_f64, 3_f64],
        ]);

        a.raw_swap(0, 2);
        assert_eq!(a, expected);
    }

    #[test]
    fn matrix_add_test() {
        const N: usize = 1000;
        const M: usize = 1000;
        let a = Matrix::new(vec![vec![1_f64; M]; N]);
        let b = Matrix::new(vec![vec![2_f64; M]; N]);
        let expected = Matrix::new(vec![vec![3_f64; M]; N]);

        let c = add(&a, &b);
        assert_eq!(c, expected);
    }

    #[test]
    fn matrix_sub_test() {
        const N: usize = 1000;
        const M: usize = 1000;
        let a = Matrix::new(vec![vec![1_f64; M]; N]);
        let b = Matrix::new(vec![vec![3_f64; M]; N]);
        let expected = Matrix::new(vec![vec![2_f64; M]; N]);

        let c = sub(&b, &a);
        assert_eq!(c, expected);
    }

    #[test]
    fn matrix_mul_test() {
        let a = Matrix::new(vec![vec![1_f64, 2_f64], vec![3_f64, 4_f64]]);
        let b = Matrix::new(vec![vec![5_f64, 6_f64], vec![7_f64, 8_f64]]);
        let expected = Matrix::new(vec![vec![19_f64, 22_f64], vec![43_f64, 50_f64]]);

        let c = mul(&a, &b);
        assert_eq!(c, expected);
    }

    #[test]
    fn gaussian_elimination_test() {
        let eps = 0.000001;
        let a = Matrix::new(vec![
            vec![2_f64, 1_f64, -1_f64],
            vec![-3_f64, -1_f64, 2_f64],
            vec![-2_f64, 1_f64, 2_f64],
        ]);
        let left = Matrix::new(vec![vec![8_f64, -11_f64, -3_f64]]);

        let expected = Matrix::new(vec![
            vec![2_f64, 1_f64, -1_f64],
            vec![0_f64, 2_f64, 1_f64],
            vec![0_f64, 0_f64, 0.25_f64],
        ]);
        let left_expected = Matrix::new(vec![vec![8_f64, 5_f64, -0.25_f64]]);

        let (a, left) = gaussian_elimination(a, left, eps);
        assert_eq!(a, expected);
        assert_eq!(left, left_expected);
    }

    #[test]
    fn system_solve_test() {
        let eps = 0.000001;
        let a = Matrix::new(vec![
            vec![2_f64, 1_f64, -1_f64],
            vec![-3_f64, -1_f64, 2_f64],
            vec![-2_f64, 1_f64, 2_f64],
        ]);
        let left = Matrix::new(vec![vec![8_f64, -11_f64, -3_f64]]);

        let expected = Matrix::new(vec![vec![2_f64, 3_f64, -1_f64]]);

        let res = system_solve(a, left, eps);
        assert_eq!(res, expected);
    }
}
