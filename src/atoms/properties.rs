use na::Matrix3xX;

use crate::atoms::new::Atoms;
use crate::constants::KB_KJPERMOLEKELVIN;

impl Atoms {
    pub fn mass_i(&self, i: usize) -> f64 {
        let type_id = self.type_ids[i] - 1;
        let mass = self.masses[type_id];
        mass
    }

    pub fn kinetic_energy(&self) -> f64 {
        let mut ek: f64 = 0.0;
        for (i, velocity) in self.velocities.column_iter().enumerate() {
            let mass = self.mass_i(i);
            ek += 0.5 * mass * velocity.norm_squared();
        }
        ek
    }

    pub fn current_temerature(&self, kinetic_energy: f64) -> f64 {
        (2.0 * kinetic_energy) / (3.0 * self.n_atoms as f64 * KB_KJPERMOLEKELVIN)
    }

    pub fn current_acceleration(&self) -> Matrix3xX<f64> {
        let mut acceleration = Matrix3xX::zeros(self.n_atoms);
        for i in 0..self.n_atoms {
            let mut a_i = acceleration.column_mut(i);
            a_i += self.forces.column(i) / self.mass_i(i);
        }
        acceleration
    }
}
