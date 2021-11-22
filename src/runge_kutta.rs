
use crate::Vec;

//  c1 | a11 a12  .. a1S
//  c2 | a21 a22  .. a2S
//  .. |  ..  ..  ..  ..
//  cS | aS1 aS2  .. aSS
// ----+-----------------
//     |  b1  b2  ..  bS
#[derive(Debug)]
pub struct ButcherTableau<const S: usize>
    where [(); S+1]: Sized
{
    a_table: [Vec<S>; S],
    b_table: Vec<S>,
    c_table: Vec<S>,
}

impl<const S: usize> ButcherTableau<S>
    where [(); S+1]: Sized
{
    // const
    fn new() -> Self {
        // let
        Self {
            a_table: [Vec::zeros(); S],
            b_table: Vec::zeros(),
            c_table: Vec::zeros(),
        }
    }
}

struct ButcherTableauBuilder<const S: usize> {
    a_nums: [Vec<S>; S],
    a_dens: [Vec<S>; S],
    b_nums: Vec<S>,
    b_dens: Vec<S>,
    c_nums: Vec<S>,
    c_dens: Vec<S>,
}

impl<const S: usize> ButcherTableauBuilder<S>
    where [(); S+1]: Sized
{
    // const
    fn build(&self) -> ButcherTableau<S> {
        let mut tableau = ButcherTableau::new();

        tableau.a_table = self.a_nums.zip(self.a_dens).map(
            |(num, den)| num.component_div(&den)
        );

        tableau.b_table = self.b_nums.component_div(&self.b_dens);
        tableau.c_table = self.c_nums.component_div(&self.c_dens);

        // dbg!(tableau)
        tableau
    }
}

impl ButcherTableau<4> {
    // const
    pub fn explicit_4() -> Self {
        let builder = ButcherTableauBuilder {
            a_nums: [
                Vec::from_iterator([0.0, 0.0, 0.0, 0.0]),
                Vec::from_iterator([1.0, 0.0, 0.0, 0.0]),
                Vec::from_iterator([0.0, 1.0, 0.0, 0.0]),
                Vec::from_iterator([0.0, 0.0, 1.0, 0.0]),
            ],
            a_dens: [
                Vec::from_iterator([1.0, 1.0, 1.0, 1.0]),
                Vec::from_iterator([2.0, 1.0, 1.0, 1.0]),
                Vec::from_iterator([1.0, 2.0, 1.0, 1.0]),
                Vec::from_iterator([1.0, 1.0, 1.0, 1.0]),
            ],
            b_nums: Vec::from_iterator([1.0, 1.0, 1.0, 1.0]),
            b_dens: Vec::from_iterator([6.0, 3.0, 3.0, 6.0]),
            c_nums: Vec::from_iterator([0.0, 1.0, 1.0, 1.0]),
            c_dens: Vec::from_iterator([1.0, 2.0, 2.0, 1.0]),
        };

        builder.build()
    }

    // const
    pub fn explicit_3_8() -> Self {
        let builder = ButcherTableauBuilder {
            a_nums: [
                Vec::from_iterator([ 0.0,  0.0, 0.0, 0.0]),
                Vec::from_iterator([ 1.0,  0.0, 0.0, 0.0]),
                Vec::from_iterator([-1.0,  1.0, 0.0, 0.0]),
                Vec::from_iterator([ 1.0, -1.0, 1.0, 0.0]),
            ],
            a_dens: [
                Vec::from_iterator([1.0, 1.0, 1.0, 1.0]),
                Vec::from_iterator([3.0, 1.0, 1.0, 1.0]),
                Vec::from_iterator([3.0, 1.0, 1.0, 1.0]),
                Vec::from_iterator([1.0, 1.0, 1.0, 1.0]),
            ],
            b_nums: Vec::from_iterator([1.0, 3.0, 3.0, 1.0]),
            b_dens: Vec::from_iterator([8.0, 8.0, 8.0, 8.0]),
            c_nums: Vec::from_iterator([0.0, 1.0, 2.0, 1.0]),
            c_dens: Vec::from_iterator([1.0, 3.0, 3.0, 1.0]),
        };

        builder.build()
    }
}

impl<const S: usize> ButcherTableau<S>
    where [(); S+1]: Sized
{
    fn k_table<const N: usize>(&self) -> [Vec<N>; S+1] {
        [Vec::zeros(); S+1]
    }

    fn sum_a_k<const N: usize>(
        &self,
        stage_num: usize,
        k_table: &[Vec<N>; S+1]
    ) -> Vec<N> {
        let mut a_k_sum = Vec::<N>::zeros();

        for idx in 0..stage_num {
            // let a = dbg!(self.a_table[stage_num]);
            // let k = dbg!(k_table[stage_num]);
            let a = self.a_table[stage_num];
            let k = k_table[stage_num];
            a_k_sum += a[idx] * k;
        }

        // dbg!(a_k_sum)
        a_k_sum
    }

    fn sum_b_k<const N: usize>(
        &self,
        k_table: &[Vec<N>; S+1],
    ) -> Vec<N> {
        let mut b_k_sum = Vec::<N>::zeros();

        for idx in 0..S {
            let b = self.b_table[idx];
            let k = k_table[idx+1];
            b_k_sum += b * k;
        }

        // dbg!(b_k_sum)
        b_k_sum
    }
}

pub fn explicit<
    const S: usize,
    const N: usize,
>(
    t_prev: &f64,
    y_prev: &Vec<N>,
    f: &dyn Fn(&f64, &Vec<N>) -> Vec<N>,
    h: f64,
    tableau: &ButcherTableau<S>
) -> (f64, Vec<N>)
    where [(); S+1]: Sized
{
    let mut k_table = tableau.k_table::<N>();

    for stage_num in 0..S {
        let del_t_s = h * tableau.c_table[stage_num];
        let del_y_s = h * tableau.sum_a_k(stage_num, &k_table);

        let t_curr_s = t_prev + del_t_s;
        let y_curr_s = y_prev + del_y_s;

        // k_table[stage_num + 1] = dbg!(f(&t_curr_s, &y_curr_s));
        k_table[stage_num + 1] = f(&t_curr_s, &y_curr_s);
    }

    let del_t = h;
    let del_y = h * tableau.sum_b_k(&k_table);

    let t_next = t_prev + del_t;
    let y_next = y_prev + del_y;

    (t_next, y_next)
}