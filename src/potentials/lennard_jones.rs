use na::Vector3;

use crate::{atoms::new::Atoms, potentials::potential::{PairPotential, PairPotentialManager, PotentialManager, Table}};

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

pub struct LennardJonesManager {
    pub table: Table,
}

impl PotentialManager for LennardJonesManager {
    fn compute_potential(&self, atoms: &mut Atoms) -> f64 {
        let mut potential_energy: f64 = 0.0;
        for i in 0..atoms.n_atoms {
            for j in (i + 1)..atoms.n_atoms {
                let mut rij = atoms.positions.column(j) - atoms.positions.column(i);

                atoms.sim_box.apply_boundary_conditions_dis(&mut rij);
                let potential = match self.get_potential_ij(atoms, i, j) {
                    Some(pot) => pot,
                    None => {
                        println!("During force calculation between {} and {} atoms, potential was missing", i + 1, j + 1);
                        continue;
                    }
                };
                if rij.norm() > potential.get_rcut() {
                    continue;
                }

                let (uij, force_ij) = potential.compute_potetial(&rij);

                potential_energy+= uij;
                {
                    let mut fi = atoms.forces.column_mut(i);
                    fi -= force_ij;

                    let mut fj = atoms.forces.column_mut(j);
                    fj += force_ij;
                }
            }
        }
        potential_energy
    }
}

impl PairPotentialManager for LennardJonesManager {
    fn with_table(table: Table) -> Self {
        Self{ table }
    }

    fn table(&self) -> &Table {
        &self.table
    }

    fn table_mut(&mut self) -> &mut Table {
        &mut self.table
    }
}
