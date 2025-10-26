use na::{DVector, Vector3};
use rand_distr::{Distribution, Normal};

use crate::atoms::new::Atoms;
use crate::constants::KB_KJPERMOLEKELVIN;

impl Atoms {
    pub fn start_velocities(&mut self, temperature: f64) {
        self.initialise_velocities(temperature);
        self.remove_drift();
        self.rescale_to_temperature(temperature);
    }

    fn initialise_velocities(&mut self, temperature: f64) {
        let mut rng = rand::rng();

        let mut sigma: f64;
        let mut normal: Normal<f64>;

        for i in 0..self.n_atoms {
            sigma = (KB_KJPERMOLEKELVIN * temperature / self.mass_i(i)).sqrt();
            normal = Normal::new(0.0, sigma).unwrap();

            self.velocities[(0, i)] = normal.sample(&mut rng);
            self.velocities[(1, i)] = normal.sample(&mut rng);
            self.velocities[(2, i)] = normal.sample(&mut rng);
        }
    }

    fn remove_drift(&mut self) {
        let mut total_mass: f64 = 0.0;
        let mut total_momentum: Vector3<f64> = Vector3::zeros();

        for i in 0..self.n_atoms {
            let a_mass = self.mass_i(i);
            total_mass += a_mass;
            total_momentum += self.velocities.column(i) * a_mass;
        }

        let velocity_cm = total_momentum / total_mass;

        let ones = DVector::from_element(self.n_atoms, 1.0);
        // (3 × 1) * (1 × n) = (3 × n)
        self.velocities -= velocity_cm * ones.transpose();
    }

    fn rescale_to_temperature(&mut self, temperature: f64) {
        let kinetic_energy = self.kinetic_energy();
        let current_temerature = self.temerature(kinetic_energy);

        let lambda = (temperature / current_temerature).sqrt();

        self.velocities *= lambda;
    }
}
