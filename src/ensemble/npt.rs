use na::{Matrix3};

use crate::{atoms::new::Atoms, constants::KB_KJPERMOLEKELVIN, math::symmetrize};

// Martyna, Tobias, Klein (1994) "Constant pressure molecular dynamics algorithms". J. Chem. Phys..
pub struct MTKBarostat {
    pub target_pressure: Matrix3<f64>,
    // barostat momentum
    pub momentum: Matrix3<f64>,
    // barostat mass
    pub w: f64,
}

impl MTKBarostat {
    pub fn new(target_pressure: Matrix3<f64>, tau: f64, n_atoms: usize, target_temp: f64) -> Self {
        let momentum = Matrix3::zeros();
        let w = ((3 * n_atoms) as f64) * KB_KJPERMOLEKELVIN * target_temp * tau.powi(2);

        Self {
            target_pressure,
            momentum,
            w,
        }
    }

    pub fn delta_momentum(&self, atoms: &Atoms, dt: f64) -> Matrix3<f64> {
        let instant_pressure = atoms.pressure_tensor();
        let delta_momentum = (instant_pressure - self.target_pressure) * (atoms.sim_box.volume() * 0.5 * dt);
        symmetrize(&delta_momentum)
    }

    pub fn scale(&self, dt: f64, velocity_scaling: bool) -> Matrix3<f64> {
        let mut eta_dot = self.momentum / self.w;
        eta_dot = symmetrize(&eta_dot);
        // let scale = mat_exp_taylor(&(-eta_dot * 0.5 * dt));
        let factor = if velocity_scaling { -0.5 } else { 1.0 };
        (eta_dot * factor * dt).exp()
    }

    pub fn kinetic_energy(&self) -> f64 {
        (self.momentum * self.momentum.transpose()).trace() / (2.0 * self.w)
    }

    pub fn potential_energy(&self, h: &Matrix3<f64>) -> f64 {
        (self.target_pressure.transpose() * h).trace()
    }
}
