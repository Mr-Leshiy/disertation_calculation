use rayon::prelude::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

#[derive(Debug, Clone, PartialEq)]
pub struct Matrix<const N: usize, const M: usize> {
    elements: Vec<Vec<f64>>,
}

impl<const N: usize, const M: usize> Matrix<N, M> {
    pub fn new(elements: Vec<Vec<f64>>) -> Self {
        assert_eq!(elements.len(), N);
        elements.par_iter().for_each(|row| assert_eq!(row.len(), M));

        Self { elements }
    }

    pub fn get_element(&self, row: usize, column: usize) -> f64 {
        self.elements[row][column]
    }

    pub fn update_element(&mut self, row: usize, column: usize, val: f64) {
        self.elements[row][column] = val;
    }

    pub fn raw_swap(&mut self, row1: usize, row2: usize) {
        if row1 != row2 {
            (0..M).into_iter().for_each(|column| {
                let tmp = self.elements[row1][column];
                self.elements[row1][column] = self.elements[row2][column];
                self.elements[row2][column] = tmp;
            });
        }
    }
}

pub fn add<const N: usize, const M: usize>(a: &Matrix<N, M>, b: &Matrix<N, M>) -> Matrix<N, M> {
    let res = (0..N)
        .into_par_iter()
        .map(|row| {
            (0..M)
                .into_par_iter()
                .map(|column| a.get_element(row, column) + b.get_element(row, column))
                .collect()
        })
        .collect();

    Matrix::<N, M>::new(res)
}

pub fn sub<const N: usize, const M: usize>(a: &Matrix<N, M>, b: &Matrix<N, M>) -> Matrix<N, M> {
    let res = (0..N)
        .into_par_iter()
        .map(|row| {
            (0..M)
                .into_par_iter()
                .map(|column| a.get_element(row, column) - b.get_element(row, column))
                .collect()
        })
        .collect();

    Matrix::<N, M>::new(res)
}

pub fn mul<const N1: usize, const M1: usize, const N2: usize, const M2: usize>(
    a: &Matrix<N1, M1>,
    b: &Matrix<N2, M2>,
) -> Matrix<N1, M2> {
    assert_eq!(M1, N2, "cannot perfom matrix multiplication");

    let res: Vec<Vec<_>> = (0..N1)
        .into_par_iter()
        .map(|row_1| {
            (0..M2)
                .into_par_iter()
                .map(|column_2| {
                    (0..M1)
                        .zip(0..N2)
                        .map(|(comunm_1, row_2)| {
                            a.get_element(row_1, comunm_1) * b.get_element(row_2, column_2)
                        })
                        .sum::<f64>()
                })
                .collect()
        })
        .collect();

    assert_eq!(res.len(), N1);
    assert_eq!(res[0].len(), M2);

    Matrix::new(res)
}

pub fn gaussian_elimination<const N: usize, const M: usize>(
    mut a: Matrix<N, M>,
    mut left: Matrix<N, 1>,
    eps: f64,
) -> (Matrix<N, M>, Matrix<N, 1>) {
    let mut h = 0;
    let mut k = 0;
    while h < N && k < M {
        let mut max = f64::abs(a.get_element(h, k));
        let mut i_max = h;
        (h..N).into_iter().for_each(|i| {
            if a.get_element(i, k) > max {
                max = a.get_element(i, k);
                i_max = i;
            }
        });
        if f64::abs(a.get_element(i_max, k)) < eps {
            k += 1;
        } else {
            a.raw_swap(h, i_max);
            left.raw_swap(h, i_max);
            (h + 1..N).into_iter().for_each(|i| {
                let f = a.get_element(i, k) / a.get_element(h, k);
                a.update_element(i, k, 0_f64);

                let left_new_val = left.get_element(i, 0) - left.get_element(h, 0) * f;
                left.update_element(i, 0, left_new_val);

                (k + 1..M).into_iter().for_each(|j| {
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

pub fn system_solve<const N: usize, const M: usize>(
    a: Matrix<N, M>,
    left: Matrix<N, 1>,
    eps: f64,
) -> Matrix<N, 1> {
    let (a, left) = gaussian_elimination(a, left, eps);
    let mut res = Matrix::<N, 1>::new(vec![vec![0_f64]; N]);

    (0..N).rev().for_each(|i| {
        let sum = (i + 1..M)
            .rev()
            .map(|j| res.get_element(j, 0) * a.get_element(i, j))
            .sum::<f64>();
        let left_val = left.get_element(i, 0) - sum;
        let val = left_val / a.get_element(i, i);
        res.update_element(i, 0, val);
    });
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matrix_raw_swap_test() {
        let mut a = Matrix::<3, 3>::new(vec![
            vec![1_f64, 2_f64, 3_f64],
            vec![4_f64, 5_f64, 6_f64],
            vec![7_f64, 8_f64, 9_f64],
        ]);
        let expected = Matrix::<3, 3>::new(vec![
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
        let a = Matrix::<N, M>::new(vec![vec![1_f64; M]; N]);
        let b = Matrix::<N, M>::new(vec![vec![2_f64; M]; N]);
        let expected = Matrix::<N, M>::new(vec![vec![3_f64; M]; N]);

        let c = add(&a, &b);
        assert_eq!(c, expected);
    }

    #[test]
    fn matrix_sub_test() {
        const N: usize = 1000;
        const M: usize = 1000;
        let a = Matrix::<N, M>::new(vec![vec![1_f64; M]; N]);
        let b = Matrix::<N, M>::new(vec![vec![3_f64; M]; N]);
        let expected = Matrix::<N, M>::new(vec![vec![2_f64; M]; N]);

        let c = sub(&b, &a);
        assert_eq!(c, expected);
    }

    #[test]
    fn matrix_mul_test() {
        let a = Matrix::<2, 2>::new(vec![vec![1_f64, 2_f64], vec![3_f64, 4_f64]]);
        let b = Matrix::<2, 2>::new(vec![vec![5_f64, 6_f64], vec![7_f64, 8_f64]]);
        let expected = Matrix::<2, 2>::new(vec![vec![19_f64, 22_f64], vec![43_f64, 50_f64]]);

        let c = mul(&a, &b);
        assert_eq!(c, expected);
    }

    #[test]
    fn gaussian_elimination_test() {
        let eps = 0.000001;
        let a = Matrix::<3, 3>::new(vec![
            vec![2_f64, 1_f64, -1_f64],
            vec![-3_f64, -1_f64, 2_f64],
            vec![-2_f64, 1_f64, 2_f64],
        ]);
        let left = Matrix::<3, 1>::new(vec![vec![8_f64], vec![-11_f64], vec![-3_f64]]);

        let expected = Matrix::<3, 3>::new(vec![
            vec![2_f64, 1_f64, -1_f64],
            vec![0_f64, 2_f64, 1_f64],
            vec![0_f64, 0_f64, 0.25_f64],
        ]);
        let left_expected = Matrix::<3, 1>::new(vec![vec![8_f64], vec![5_f64], vec![-0.25_f64]]);

        let (a, left) = gaussian_elimination(a, left, eps);
        assert_eq!(a, expected);
        assert_eq!(left, left_expected);
    }

    #[test]
    fn system_solve_test() {
        let eps = 0.000001;
        let a = Matrix::<3, 3>::new(vec![
            vec![2_f64, 1_f64, -1_f64],
            vec![-3_f64, -1_f64, 2_f64],
            vec![-2_f64, 1_f64, 2_f64],
        ]);
        let left = Matrix::<3, 1>::new(vec![vec![8_f64], vec![-11_f64], vec![-3_f64]]);

        let expected = Matrix::<3, 1>::new(vec![vec![2_f64], vec![3_f64], vec![-1_f64]]);

        let res = system_solve(a, left, eps);
        assert_eq!(res, expected);
    }
}
