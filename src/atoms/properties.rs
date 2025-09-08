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

    pub fn current_temerature(&self) -> f64 {
        let kinetic_energy = self.kinetic_energy();
        (2.0 * kinetic_energy ) / (3.0 * self.n_atoms as f64 * KB_KJPERMOLEKELVIN)
    }
}