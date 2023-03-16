use rayon::prelude::{IntoParallelIterator, ParallelIterator};

pub trait MatrixTrait<const N: usize, const M: usize>: Sync {
    fn get_element_impl(&self, i: usize, j: usize) -> f64;

    fn get_element(&self, i: usize, j: usize) -> f64 {
        assert!(i < N);
        assert!(j < M);

        self.get_element_impl(i, j)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Matrix<const N: usize, const M: usize> {
    elements: Vec<Vec<f64>>,
}

impl<const N: usize, const M: usize> Matrix<N, M> {
    pub fn new(elements: Vec<Vec<f64>>) -> Self {
        assert_eq!(elements.len(), N);
        assert_eq!(elements[0].len(), M);

        Self { elements }
    }
}

impl<const N: usize, const M: usize> MatrixTrait<N, M> for Matrix<N, M> {
    fn get_element_impl(&self, i: usize, j: usize) -> f64 {
        self.elements[i][j]
    }
}

pub fn add<const N: usize, const M: usize, T1: MatrixTrait<N, M>, T2: MatrixTrait<N, M>>(
    a: &T1,
    b: &T2,
) -> Matrix<N, M> {
    let res = (0..N)
        .into_par_iter()
        .map(|i| {
            (0..M)
                .into_par_iter()
                .map(move |j| a.get_element(i, j) + b.get_element(i, j))
                .collect()
        })
        .collect();

    Matrix::<N, M>::new(res)
}

pub fn sub<const N: usize, const M: usize, T1: MatrixTrait<N, M>, T2: MatrixTrait<N, M>>(
    a: &T1,
    b: &T2,
) -> Matrix<N, M> {
    let res = (0..N)
        .into_par_iter()
        .map(|i| {
            (0..M)
                .into_par_iter()
                .map(move |j| a.get_element(i, j) - b.get_element(i, j))
                .collect()
        })
        .collect();

    Matrix::<N, M>::new(res)
}

// pub fn mul<
//     const N1: usize,
//     const M1: usize,
//     T1: MatrixTrait<N1, M1>,
//     const N2: usize,
//     const M2: usize,
//     T2: MatrixTrait<N2, M2>,
// >(
//     a: &T1,
//     b: &T2,
// ) -> Matrix<N1, M2> {
//     assert_eq!(M1, N2, "cannot perfom matrix multiplication");

//     let res = Matrix::<N1, M2>::new([[0_f64; M2]; N1]);
//     res
// }

#[cfg(test)]
mod tests {
    use super::*;

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
}
