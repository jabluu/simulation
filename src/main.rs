#![feature(
    const_for,
    const_fn_floating_point_arithmetic,
    const_mut_refs,
    generic_const_exprs,
    array_zip,
)]
// #![allow(unused)]

use nalgebra as na;
mod runge_kutta;

type Vec<const C: usize> = na::SVector<f64, C>;
type Mat<const R: usize, const C: usize> = na::SMatrix<f64, R, C>;

type Vec1 = Vec<1>;
type Vec2 = Vec<2>;
type Vec3 = Vec<3>;
type Vec4 = Vec<4>;
type Vec5 = Vec<5>;
type Vec6 = Vec<6>;

type Mat3 = Mat<3,3>;

#[derive(Debug)]
struct RigidBody {
    lin_mass: f64,
    rot_mass: Mat3,
    lin_pos: Vec3,
    lin_vel: Vec3,
    rot_pos: Vec4,
    rot_vel: Vec3,
}

impl std::fmt::Display for RigidBody {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{} {} {} {} {} {} {} {} {} {} {} {}",
            self.lin_pos[0], self.lin_pos[1], self.lin_pos[2],
            self.lin_vel[0], self.lin_vel[1], self.lin_vel[2],
            self.rot_pos[0], self.rot_pos[1], self.rot_pos[2],
            self.rot_vel[0], self.rot_vel[1], self.rot_vel[2],
        )
    }
}

const GRAVITATIONAL_CONSTANT: f64 = 6.67430e-11; // m^3 kg^-1 s^-2

fn newtonian_gravity(
    rigidbody1: &RigidBody,
    rigidbody2: &RigidBody
) -> Vec3 {
    let G = GRAVITATIONAL_CONSTANT;

    let m1 = rigidbody1.lin_mass;
    let m2 = rigidbody2.lin_mass;

    let r21 = rigidbody2.lin_pos - rigidbody1.lin_pos;
    let r21_unit = r21.normalize();
    let r21_norm_sq = r21.norm_squared();

    let f21 = -G * m1 * m2 / r21_norm_sq * r21_unit;

    f21
}

fn inverse_symmetric(m: &Mat3) -> Mat3 {
    let den =
        m[(0,0)] * m[(1,1)] * m[(2,2)]
        + 2.0 * m[(0,1)] * m[(0,2)] * m[(1,2)]
        - m[(0,0)] * m[(1,2)].powi(2)
        - m[(1,1)] * m[(0,2)].powi(2)
        - m[(2,2)] * m[(0,1)].powi(2);

    let m00 = m[(1,1)] * m[(2,2)] - m[(1,2)].powi(2);
    let m11 = m[(0,0)] * m[(2,2)] - m[(0,2)].powi(2);
    let m22 = m[(0,0)] * m[(1,1)] - m[(0,1)].powi(2);

    let m01 = m[(0,2)] * m[(1,2)] - m[(0,1)] * m[(2,2)];
    let m02 = m[(0,1)] * m[(1,2)] - m[(0,2)] * m[(1,1)];
    let m12 = m[(0,1)] * m[(0,2)] - m[(0,0)] * m[(1,2)];

    let inv_m = Mat::from_iterator([
        m00, m01, m02,
        m01, m11, m12,
        m02, m12, m22,
    ]) / den;

    inv_m
}

fn inverse_principal(m: &Mat3) -> Mat3 {
    Mat::from_iterator([
        1.0 / m[(0,0)], 0.0, 0.0,
        0.0, 1.0 / m[(1,1)], 0.0,
        0.0, 0.0, 1.0 / m[(2,2)],
    ])
}

fn newton_euler_dynamics(
    rigidbody: &RigidBody,
    external_forces: &Vec3,
    external_moments: &Vec3,
) -> (Vec3, Vec3) {
    let RigidBody {
        lin_mass,
        rot_mass,
        rot_vel,
        ..
    } = rigidbody;

    let inv_lin_mass = 1.0 / lin_mass;
    let inv_rot_mass = inverse_symmetric(&rot_mass);

    let rot_mass_vel = rot_mass * rot_vel;

    let lin_mom_rate = external_forces;
    let rot_mom_rate = external_moments - rot_vel.cross(&rot_mass_vel);

    let lin_acc = inv_lin_mass * lin_mom_rate;
    let rot_acc = inv_rot_mass * rot_mom_rate;

    (lin_acc, rot_acc)
}

fn main() {
    let tableau = runge_kutta::ButcherTableau::explicit_4();

    let mut time = 0.0; // [s]
    let rate = 100.0; // [Hz]

    let mut rigidbody1 = RigidBody {
        lin_mass: 1.0,
        rot_mass: Mat::identity(),
        lin_pos: Vec::from_iterator([-1.0, 0.0, 0.0]),
        lin_vel: Vec::from_iterator([0.0, 0.0, 0.0]),
        rot_pos: Vec::zeros(),
        rot_vel: Vec::zeros(),
    };

    let mut rigidbody2 = RigidBody {
        lin_mass: 1.0,
        rot_mass: Mat::identity(),
        lin_pos: Vec::from_iterator([1.0, 0.0, 0.0]),
        lin_vel: Vec::from_iterator([0.0, 0.0, 0.0]),
        rot_pos: Vec::zeros(),
        rot_vel: Vec::zeros(),
    };

    loop {
        let f21 = newtonian_gravity(&rigidbody1, &rigidbody2);

        let (lin_acc, rot_acc) = newton_euler_dynamics(
            &rigidbody1,
            &f21,
            &Vec::zeros(),
        );

        let (new_time, new_lin_pos) = runge_kutta::explicit(
            &time,
            &rigidbody1.lin_pos,
            &|t, y| rigidbody1.lin_vel,
            1.0 / rate,
            &tableau
        );

        let (new_time, new_lin_vel) = runge_kutta::explicit(
            &time,
            &rigidbody1.lin_vel,
            &|t, y| lin_acc,
            1.0 / rate,
            &tableau
        );

        rigidbody1.lin_pos = new_lin_pos;
        rigidbody1.lin_vel = new_lin_vel;

        let (lin_acc, rot_acc) = newton_euler_dynamics(
            &rigidbody2,
            &(-f21),
            &Vec::zeros(),
        );

        let (new_time, new_lin_pos) = runge_kutta::explicit(
            &time,
            &rigidbody2.lin_pos,
            &|t, y| rigidbody2.lin_vel,
            1.0 / rate,
            &tableau
        );

        let (new_time, new_lin_vel) = runge_kutta::explicit(
            &time,
            &rigidbody2.lin_vel,
            &|t, y| lin_acc,
            1.0 / rate,
            &tableau
        );

        rigidbody2.lin_pos = new_lin_pos;
        rigidbody2.lin_vel = new_lin_vel;

        time = new_time;

        println!("{} {} {}", time, rigidbody1, rigidbody2);

        if time > 10.0 {
            break
        }

        // let gravitational_acceleration = Vec3::from_iterator([
        //     0.0,
        //     0.0,
        //     0.0,
        // ]);

        // let acceleration = Vec3::from_iterator([
        //     gravitational_acceleration[0],
        //     gravitational_acceleration[1],
        //     gravitational_acceleration[2],
        // ]);
    }
}
