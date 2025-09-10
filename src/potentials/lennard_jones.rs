use na::Vector3;

use crate::potentials::potential::PairPotential;

pub struct LennardJones {
    epsilon: f64,
    sigma: f64,
    rcut: f64,
    shift: bool,
}

impl LennardJones {
    pub fn new(epsilon: f64, sigma: f64, rcut: f64, shift: bool) -> Self {
        Self {
            epsilon,
            sigma,
            rcut,
            shift
        }
    }
}

impl PairPotential for LennardJones {
    fn compute_potetial(&self, rij: &Vector3<f64>) -> (f64, Vector3<f64>) {
        let rij2 = rij.norm_squared();
        let inv_rij2 = 1.0 / rij2;
        let vanderwaals_attraction = (self.sigma.powi(2) * inv_rij2).powi(3);
        let lj_repulsion = vanderwaals_attraction.powi(2);

        let mut potential_energy = 4.0 * self.epsilon * (lj_repulsion - vanderwaals_attraction);

        let force = 24.0 * self.epsilon * (2.0 * lj_repulsion - vanderwaals_attraction) * inv_rij2 * rij;

        if self.shift {
            let cutoff_inv2 = (self.sigma / self.rcut).powi(2);
            let cutoff_attraction = cutoff_inv2.powi(3);
            let cutoff_repulsion = cutoff_attraction.powi(2);

            let u_cutoff = 4.0 * self.epsilon * (cutoff_repulsion - cutoff_attraction);

            potential_energy -= u_cutoff;
        }

        (potential_energy, force)
    }

    fn get_rcut(&self) -> f64 {
        self.rcut
    }
}
