use na::Matrix3xX;

use crate::atoms::new::Atoms;
use crate::constants::KB_KJPERMOLEKELVIN;


impl Atoms {
    pub fn kinetic_energy(&self) -> f64 {
        self.velocities
            .column_iter()
            .zip(self.masses.iter())
            .map(|(velocity, &mass)| 0.5 * mass * velocity.norm_squared())
            .sum()
    }

    pub fn current_temerature(&self, kinetic_energy: f64) -> f64 {
        (2.0 * kinetic_energy ) / (3.0 * self.n_atoms as f64 * KB_KJPERMOLEKELVIN)
    }

    pub fn current_acceleration(&self) -> Matrix3xX<f64> {
        let mut acceleration = Matrix3xX::zeros(self.n_atoms);
        for i in 0..self.n_atoms {
            let mut a_i = acceleration.column_mut(i);
            a_i += self.forces.column(i) / self.masses[i]; 
        }
        acceleration
    }
}
