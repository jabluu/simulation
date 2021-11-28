
use crate::Vec;

//  c1 | a11 a12  .. a1S
//  c2 | a21 a22  .. a2S
//  .. |  ..  ..  ..  ..
//  cS | aS1 aS2  .. aSS
// ----+-----------------
//     |  b1  b2  ..  bS
#[derive(Debug)]
pub struct ButcherTableau<const S: usize> {
    a_table: [[f64; S]; S], // [Vec<S>; S],
    b_table: [f64; S], // Vec<S>,
    c_table: [f64; S], // Vec<S>,
}

impl<const S: usize> ButcherTableau<S> {
    // const
    fn new() -> Self {
        // let
        Self {
            a_table: [[0.0; S]; S], // [[Vec::zeros()]; S],
            b_table: [0.0; S], // Vec::zeros(),
            c_table: [0.0; S], // Vec::zeros(),
        }
    }
}

impl ButcherTableau<4> {
    // const
    pub fn explicit4() -> Self {
        ButcherTableau {
            a_table: [
                [0.0 / 1.0, 0.0 / 1.0, 0.0 / 1.0, 0.0 / 1.0],
                [1.0 / 2.0, 0.0 / 1.0, 0.0 / 1.0, 0.0 / 1.0],
                [0.0 / 1.0, 1.0 / 2.0, 0.0 / 1.0, 0.0 / 1.0],
                [0.0 / 1.0, 0.0 / 1.0, 1.0 / 1.0, 0.0 / 1.0],
            ],
            b_table: [1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0],
            c_table: [0.0 / 1.0, 1.0 / 2.0, 1.0 / 2.0, 1.0 / 1.0],
        }
    }

    // const
    pub fn explicit38() -> Self {
        ButcherTableau {
            a_table: [
                [ 0.0 / 1.0,  0.0 / 1.0, 0.0 / 1.0, 0.0 / 1.0],
                [ 1.0 / 3.0,  0.0 / 1.0, 0.0 / 1.0, 0.0 / 1.0],
                [-1.0 / 3.0,  1.0 / 1.0, 0.0 / 1.0, 0.0 / 1.0],
                [ 1.0 / 1.0, -1.0 / 1.0, 1.0 / 1.0, 0.0 / 1.0],
            ],
            b_table: [1.0 / 8.0, 3.0 / 8.0, 3.0 / 8.0, 1.0 / 8.0],
            c_table: [0.0 / 1.0, 1.0 / 3.0, 2.0 / 3.0, 1.0 / 1.0],
        }
    }
}

#[derive(Debug)]
pub struct RungeKutta<const N: usize, const S: usize> {
    tableau: ButcherTableau<S>,
    k_table: [Vec<N>; S],
}

impl<const N: usize> RungeKutta<N, 4> {
    pub fn explicit4() -> Self {
        RungeKutta {
            k_table: [Vec::zeros(); 4],
            tableau: ButcherTableau::explicit4(),
        }
    }

    pub fn explicit38() -> Self {
        RungeKutta {
            k_table: [Vec::zeros(); 4],
            tableau: ButcherTableau::explicit38(),
        }
    }
}


impl<const N: usize, const S: usize> RungeKutta<N, S> {
    pub fn integrate(
        &self,
        t: f64,
        h: f64,
        y: &Vec<N>,
        f: &dyn Fn(f64, Vec<N>) -> Vec<N>,
    ) -> Vec<N> {
        let RungeKutta {
            tableau: ButcherTableau {
                a_table,
                b_table,
                c_table,
                ..
            },
            mut k_table,
            ..
        } = self;

        let mut sum_b_k = Vec::zeros();

        for stage_num in 0..S {
            let mut sum_a_k = Vec::zeros();

            for idx in 0..stage_num {
                sum_a_k += a_table[stage_num][idx] * k_table[idx];
            }

            k_table[stage_num] = f(
                t + h * c_table[stage_num],
                y + h * sum_a_k,
            );

            sum_b_k += b_table[stage_num] * k_table[stage_num];
        }

        y + h * sum_b_k
    }
}