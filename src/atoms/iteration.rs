use crate::atoms::new::Atoms;
use crate::potentials::lennard_jones::Potentials;

impl Atoms {
    pub fn iterate(&mut self, sigma: f64, epsilon: f64, rcut: f64, shift: bool) -> f64 {
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
}
