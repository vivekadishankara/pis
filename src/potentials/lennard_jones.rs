use std::collections::HashSet;

use na::Vector3;

use crate::{atoms::new::Atoms, potentials::potential::{PairPotential, PairPotentialManager, PotentialManager, Table}};
use crate::atoms::neighbour_list::FORWARD_NEIGHBOUR_OFFSETS;

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
    fn compute_potential(&self, rij: &Vector3<f64>) -> (f64, Vector3<f64>) {
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
                if i == 0 {
                    print!("{} {}, ", i, j); // =========================================================
                }

                let (uij, force_ij) = potential.compute_potential(&rij);

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

pub struct LennardJonesVerletManager {
    pub table: Table,
}

impl PotentialManager for LennardJonesVerletManager {
    fn compute_potential (&self, atoms: &mut Atoms) -> f64 {
        let max_rcut = self.max_rcut();

        let (nx, ny, nz) = atoms.divide_into_cells(max_rcut);
        let shift = [-1, 0, 1];
        let cells = atoms.rcut_cells(nx, ny, nz);

        let mut neighbour_cells: HashSet<usize> = HashSet::with_capacity(27);
        let mut potential_energy: f64 = 0.0;

        for cx_i in 0..nx {
            for cy_i in 0..ny {
                for cz_i in 0..nz {
                    let current_cell = Atoms::cell_index(cx_i, cy_i, cz_i, nx, ny);
                    println!("current cell {current_cell}"); // =========================================================
                    let current_cell_atoms = &cells[current_cell];

                    for dx in &shift {
                        for dy in &shift {
                            for dz in &shift {
                                let cx_j = (cx_i as isize + dx).rem_euclid(nx as isize) as usize;
                                let cy_j = (cy_i as isize + dy).rem_euclid(ny as isize) as usize;
                                let cz_j = (cz_i as isize + dz).rem_euclid(nz as isize) as usize;

                                let a_neighbour_cell = Atoms::cell_index(cx_j, cy_j, cz_j, nx, ny);
                                if neighbour_cells.contains(&a_neighbour_cell) { continue; }
                                neighbour_cells.insert(a_neighbour_cell);
                                let neighbour_cell_atoms = &cells[a_neighbour_cell];

                                for &i in current_cell_atoms.iter() {
                                    for &j in neighbour_cell_atoms.iter() {
                                        if i <= j { continue; }
                                        let mut rij = atoms.positions.column(j) - atoms.positions.column(i);
                                        atoms.sim_box.apply_boundary_conditions_dis(&mut rij);

                                        let potential = match self.get_potential_ij(atoms, i,j) {
                                            Some(pot) => pot,
                                            None => {
                                                println!("During force calculation between {} and {} atoms, potential was missing", i + 1, j + 1);
                                                continue;
                                            }                                            
                                        };

                                        if rij.norm() > potential.get_rcut() { continue; }
                                        let (uij, force_ij) = potential.compute_potential(&rij);
                                        potential_energy += uij;
                                        {
                                            let mut fi = atoms.forces.column_mut(i);
                                            fi -= force_ij;

                                            let mut fj = atoms.forces.column_mut(j);
                                            fj += force_ij;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    neighbour_cells.clear();
                }
            }
        }
        potential_energy
    }
}

impl PairPotentialManager for LennardJonesVerletManager {
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

pub struct LennardJonesVerletOffsetManager {
    pub table: Table,
}

impl PotentialManager for LennardJonesVerletOffsetManager {
    fn compute_potential(&self, atoms: &mut Atoms) -> f64 {
        let max_rcut = self.max_rcut();
        let (nx, ny, nz) = atoms.divide_into_cells(max_rcut);
        let cells = atoms.rcut_cells(nx, ny, nz);

        let mut neighbour_cells: HashSet<usize> = HashSet::with_capacity(27);
        let mut potential_energy: f64 = 0.0;

        for cx_i in 0..nx {
            for cy_i in 0..ny {
                for cz_i in 0..nz {
                    let current_cell = Atoms::cell_index(cx_i, cy_i, cz_i, nx, ny);
                    let current_cell_atoms = &cells[current_cell];

                    for offset in FORWARD_NEIGHBOUR_OFFSETS.iter() {
                        let cx_j = (cx_i as isize + offset.dx).rem_euclid(nx as isize) as usize;
                        let cy_j = (cy_i as isize + offset.dy).rem_euclid(ny as isize) as usize;
                        let cz_j = (cz_i as isize + offset.dz).rem_euclid(nz as isize) as usize;

                        let a_neighbour_cell = Atoms::cell_index(cx_j, cy_j, cz_j, nx, ny);
                        if neighbour_cells.contains(&a_neighbour_cell) { continue; }
                        neighbour_cells.insert(a_neighbour_cell);

                        let neighbour_cell_atoms = &cells[a_neighbour_cell];

                        for &i in current_cell_atoms.iter() {
                            for &j in neighbour_cell_atoms.iter() {
                                if current_cell == a_neighbour_cell && i<=j { continue; }
                                let mut rij = atoms.positions.column(j) - atoms.positions.column(i);
                                atoms.sim_box.apply_boundary_conditions_dis(&mut rij);

                                let potential = match self.get_potential_ij(atoms, i,j) {
                                    Some(pot) => pot,
                                    None => {
                                        println!("During force calculation between {} and {} atoms, potential was missing", i + 1, j + 1);
                                        continue;
                                    }                                            
                                };

                                if rij.norm() > potential.get_rcut() { continue; }
                                let (uij, force_ij) = potential.compute_potential(&rij);
                                potential_energy += uij;
                                {
                                    let mut fi = atoms.forces.column_mut(i);
                                    fi -= force_ij;

                                    let mut fj = atoms.forces.column_mut(j);
                                    fj += force_ij;
                                }
                            }
                        }
                    }
                    neighbour_cells.clear();
                }
            }
        }
        potential_energy
    }
}

impl PairPotentialManager for LennardJonesVerletOffsetManager {
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
