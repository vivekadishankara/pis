use na::{DVector, Matrix3xX};

use crate::atoms::new::Atoms;
use crate::potentials::lennard_jones::Potentials;

impl Atoms {
    pub fn compute_potential_and_forces(
        &mut self, 
        sigma: f64, 
        epsilon: f64, 
        rcut: f64, 
        shift: bool
    ) -> f64 {
        let mut potential_energy: f64 = 0.0;
        for i in 0..self.n_atoms {
            for j in (i + 1)..self.n_atoms {
                let rij = self.positions.column(j) - self.positions.column(i);

                let rij = self.sim_box.apply_boundary_conditions(&rij);

                if rij.norm() > rcut {
                    continue;
                }

                let (uij, force_ij) = Potentials::lennard_jones(&rij, sigma, epsilon, rcut, shift);

                potential_energy+= uij;
                {
                    let mut fi = self.forces.column_mut(i);
                    fi -= force_ij;

                    let mut fj = self.forces.column_mut(j);
                    fj += force_ij;
                }
            }
        }

        potential_energy
    }

    pub fn verlet_step(
        &mut self, 
        dt: f64, 
        sigma: f64, 
        epsilon: f64, 
        rcut: f64, 
        shift: bool
    ) -> f64 {

        let a_t = self.current_acceleration();

        self.positions += &self.velocities * dt + &a_t * 0.5 * dt.powi(2);

        for mut r_i in self.positions.column_iter_mut() {
            let wrapped_r_i = self.sim_box.apply_boundary_conditions(&r_i.clone_owned());
            r_i.copy_from(&wrapped_r_i);
        }
        
        self.forces = Matrix3xX::zeros(self.n_atoms);

        let potential_energy = self.compute_potential_and_forces(sigma, epsilon, rcut, shift);

        let a_tdt = self.current_acceleration();

        self.velocities += (a_t + a_tdt) * 0.5 * dt;

        potential_energy

    }

    pub fn run(
        &mut self, 
        dt: f64, 
        sigma: f64, 
        epsilon: f64, 
        rcut: f64, 
        shift: bool, 
        time_steps: usize
    ) -> DVector<f64> {
        let mut potential_energies: DVector<f64> = DVector::zeros(time_steps + 1);
        let first_potential = self.compute_potential_and_forces(sigma, epsilon, rcut, shift);
        potential_energies[0] = first_potential;
        for i in 0..time_steps {
            let step_potential = self.verlet_step(dt, sigma, epsilon, rcut, shift);
            potential_energies[i + 1] = step_potential;
        }
        potential_energies
    }
}
