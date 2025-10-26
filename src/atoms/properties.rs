use na::{Matrix3, Matrix3xX};

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

    pub fn temerature(&self, kinetic_energy: f64) -> f64 {
        (2.0 * kinetic_energy) / (self.degress_of_freedom() as f64 * KB_KJPERMOLEKELVIN)
    }

    pub fn current_acceleration(&self) -> Matrix3xX<f64> {
        let mut acceleration = Matrix3xX::zeros(self.n_atoms);
        for i in 0..self.n_atoms {
            let mut a_i = acceleration.column_mut(i);
            a_i += self.forces.column(i) / self.mass_i(i);
        }
        acceleration
    }

    pub fn degress_of_freedom(&self) -> usize {
        3 * self.n_atoms
    }

    pub fn kinetic_tensor(&self) -> Matrix3<f64> {
        let mut kinetic_tensor: Matrix3<f64> = Matrix3::zeros();
        for i in 0..self.n_atoms {
            let mass = self.mass_i(i);
            kinetic_tensor += mass * self.velocities.column(i) * self.velocities.column(i).transpose();
        }
        kinetic_tensor
    }

    pub fn virial_tensor(&self) -> Matrix3<f64> {
        &self.positions * self.forces.transpose()
    }

    pub fn pressure_tensor(&self) -> Matrix3<f64> {
        let kinetic_tensor = self.kinetic_tensor();
        let virial_tensor = self.virial_tensor();
        let volume = self.sim_box.volume();

        (kinetic_tensor + virial_tensor) / volume
    }

    pub fn pressure(&self, kinetic_energy: f64) -> f64 {
        let virial = self.virial_tensor().trace();
        let volume = self.sim_box.volume();
        (2.0 * kinetic_energy + virial) / (3.0 * volume)
    }
}
