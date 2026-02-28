#![allow(dead_code)]
use std::{collections::HashSet, ops::AddAssign, sync::Arc};

use na::{DVector, Matrix3xX, Vector3};
use parking_lot::RwLock;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::atoms::neighbour_list::FORWARD_NEIGHBOUR_OFFSETS;
use crate::{
    atoms::new::Atoms,
    potentials::potential::{PairPotential, PairPotentialManager, PotentialManager, Table},
};

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
            shift,
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

        let force =
            24.0 * self.epsilon * (2.0 * lj_repulsion - vanderwaals_attraction) * inv_rij2 * rij;

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

pub struct LJManager {
    pub table: Table,
}

impl PotentialManager for LJManager {
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
        potential_energy
    }
}

impl_pair_potential_manager!(LJManager);

pub struct LJVerletManager {
    pub table: Table,
}

impl PotentialManager for LJVerletManager {
    fn compute_potential(&self, atoms: &mut Atoms) -> f64 {
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
                    let current_cell_atoms = &cells[current_cell];

                    for dx in &shift {
                        for dy in &shift {
                            for dz in &shift {
                                let cx_j = (cx_i as isize + dx).rem_euclid(nx as isize) as usize;
                                let cy_j = (cy_i as isize + dy).rem_euclid(ny as isize) as usize;
                                let cz_j = (cz_i as isize + dz).rem_euclid(nz as isize) as usize;

                                let a_neighbour_cell = Atoms::cell_index(cx_j, cy_j, cz_j, nx, ny);
                                if neighbour_cells.contains(&a_neighbour_cell) {
                                    continue;
                                }
                                neighbour_cells.insert(a_neighbour_cell);
                                let neighbour_cell_atoms = &cells[a_neighbour_cell];

                                for &i in current_cell_atoms.iter() {
                                    for &j in neighbour_cell_atoms.iter() {
                                        if i <= j {
                                            continue;
                                        }
                                        let mut rij =
                                            atoms.positions.column(j) - atoms.positions.column(i);
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

impl_pair_potential_manager!(LJVerletManager);

pub struct LJVOffsetManager {
    pub table: Table,
}

impl PotentialManager for LJVOffsetManager {
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
                        if neighbour_cells.contains(&a_neighbour_cell) {
                            continue;
                        }
                        neighbour_cells.insert(a_neighbour_cell);

                        let neighbour_cell_atoms = &cells[a_neighbour_cell];

                        for &i in current_cell_atoms.iter() {
                            for &j in neighbour_cell_atoms.iter() {
                                if current_cell == a_neighbour_cell && i <= j {
                                    continue;
                                }
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

impl_pair_potential_manager!(LJVOffsetManager);

pub struct LJVParallelManager {
    pub table: Table,
}

impl PotentialManager for LJVParallelManager {
    fn compute_potential(&self, atoms: &mut Atoms) -> f64 {
        let max_rcut = self.max_rcut();
        let (nx, ny, nz) = atoms.divide_into_cells(max_rcut);
        let cells = atoms.rcut_cells(nx, ny, nz);

        let cell_indices = Atoms::cell_indices_for_parallel(nx, ny, nz);

        let n_threads = rayon::current_num_threads();

        let thread_forces: Vec<Matrix3xX<f64>> = (0..n_threads)
            .map(|_| Matrix3xX::zeros(atoms.n_atoms))
            .collect();
        let thread_forces = std::sync::Arc::new(parking_lot::RwLock::new(thread_forces));

        let thread_potential_energy: DVector<f64> = DVector::zeros(n_threads);
        let thread_potential_energy =
            std::sync::Arc::new(parking_lot::RwLock::new(thread_potential_energy));

        let thread_neighbour_cells: Vec<HashSet<usize>> =
            (0..n_threads).map(|_| HashSet::new()).collect();
        let thread_neighbour_cells =
            std::sync::Arc::new(parking_lot::RwLock::new(thread_neighbour_cells));

        cell_indices.par_iter().for_each(|&(cx_i, cy_i, cz_i)| {
            let tid = rayon::current_thread_index().unwrap();
            let mut tf = thread_forces.write();
            let mut tpe = thread_potential_energy.write();
            let mut tnc = thread_neighbour_cells.write();

            let current_cell = Atoms::cell_index(cx_i, cy_i, cz_i, nx, ny);
            let current_cell_atoms = &cells[current_cell];

            for offset in FORWARD_NEIGHBOUR_OFFSETS.iter() {
                let cx_j = (cx_i as isize + offset.dx).rem_euclid(nx as isize) as usize;
                let cy_j = (cy_i as isize + offset.dy).rem_euclid(ny as isize) as usize;
                let cz_j = (cz_i as isize + offset.dz).rem_euclid(nz as isize) as usize;

                let a_neighbour_cell = Atoms::cell_index(cx_j, cy_j, cz_j, nx, ny);
                if tnc[tid].contains(&a_neighbour_cell) { continue; }
                tnc[tid].insert(a_neighbour_cell);
                let neighbour_cell_atoms = &cells[a_neighbour_cell];

                for &i in current_cell_atoms.iter() {
                    for &j in neighbour_cell_atoms.iter() {
                        if current_cell == a_neighbour_cell && i<=j { continue; }
                        let mut rij = atoms.positions.column(j) - atoms.positions.column(i);
                        atoms.sim_box.apply_boundary_conditions_dis(&mut rij);

                        let potential = match self.get_potential_ij(atoms, i, j) {
                            Some(pot) => pot,
                            None => {
                                println!("During force calculation between {} and {} atoms, potential was missing", i + 1, j + 1);
                                continue;
                            }
                        };

                        if rij.norm() > potential.get_rcut() { continue; }

                        let (uij, force_ij) = potential.compute_potential(&rij);

                        tpe[tid] += uij;
                        {
                            let mut fi = tf[tid].column_mut(i);
                            fi -= force_ij;

                            let mut fj = tf[tid].column_mut(j);
                            fj += force_ij;
                        }
                    }
                }
            }
            tnc[tid].clear();
        });

        let tf = thread_forces.read();

        for i in 0..n_threads {
            atoms.forces += &tf[i];
        }

        let tpe = thread_potential_energy.read();
        tpe.sum()
    }
}

impl_pair_potential_manager!(LJVParallelManager);

pub struct LJVPBuildListManager {
    pub table: Table,
}

impl LJVPBuildListManager {
    pub fn build_neighbour_list(&self, atoms: &mut Atoms) -> Vec<Vec<usize>> {
        let max_rcut = self.max_rcut();
        let (nx, ny, nz) = atoms.divide_into_cells(max_rcut);
        let total_cells = nx * ny * nz;
        let shift = [-1, 0, 1];

        let cells = atoms.rcut_cells(nx, ny, nz);

        let cell_indices = Atoms::cell_indices_for_parallel(nx, ny, nz);

        let n_threads = rayon::current_num_threads();
        let neighbour_list_len = (2 * atoms.n_atoms) / total_cells;
        let thread_neighbour_list: Vec<Vec<usize>> = (0..atoms.n_atoms)
            .map(|_| Vec::with_capacity(neighbour_list_len))
            .collect();
        let thread_neighbour_list =
            std::sync::Arc::new(parking_lot::RwLock::new(thread_neighbour_list));

        let thread_neighbour_cells: Vec<HashSet<usize>> =
            (0..n_threads).map(|_| HashSet::new()).collect();
        let thread_neighbour_cells =
            std::sync::Arc::new(parking_lot::RwLock::new(thread_neighbour_cells));

        cell_indices.par_iter().for_each(|&(cx_i, cy_i, cz_i)| {
            let tid = rayon::current_thread_index().unwrap();

            let mut tnc = thread_neighbour_cells.write();
            let mut tnl = thread_neighbour_list.write();

            let current_cell = Atoms::cell_index(cx_i, cy_i, cz_i, nx, ny);
            let current_cell_atoms = &cells[current_cell];

            for &dx in &shift {
                for &dy in &shift {
                    for &dz in &shift {
                        let cx_j = (cx_i as isize + dx).rem_euclid(nx as isize) as usize;
                        let cy_j = (cy_i as isize + dy).rem_euclid(ny as isize) as usize;
                        let cz_j = (cz_i as isize + dz).rem_euclid(nz as isize) as usize;

                        let a_neighbour_cell = Atoms::cell_index(cx_j, cy_j, cz_j, nx, ny);

                        if tnc[tid].contains(&a_neighbour_cell) { continue; }
                        tnc[tid].insert(a_neighbour_cell);
                        let neighbour_cell_atoms = &cells[a_neighbour_cell];

                        for &i in current_cell_atoms.iter() {
                            for &j in neighbour_cell_atoms.iter() {
                                if current_cell == a_neighbour_cell && i==j { continue; }
                                let mut rij = atoms.positions.column(j) - atoms.positions.column(i);
                                atoms.sim_box.apply_boundary_conditions_dis(&mut rij);

                                let potential = match self.get_potential_ij(atoms, i, j) {
                                    Some(pot) => pot,
                                    None => {
                                        println!("During force calculation between {} and {} atoms, potential was missing", i + 1, j + 1);
                                        continue;
                                    }
                                };

                                if rij.norm() > potential.get_rcut() { continue; }
                                tnl[i].push(j);
                            }
                        }
                    }
                }
                tnc[tid].clear();
            }
        });
        let tnl = thread_neighbour_list.read();
        tnl.to_vec()
    }
}

impl PotentialManager for LJVPBuildListManager {
    fn compute_potential(&self, atoms: &mut Atoms) -> f64 {
        let neighbour_list = self.build_neighbour_list(atoms);

        let thread_forces: Matrix3xX<f64> = Matrix3xX::zeros(atoms.n_atoms);
        let thread_forces = Arc::new(RwLock::new(thread_forces));

        let thread_potential_energy: DVector<f64> = DVector::zeros(atoms.n_atoms);
        let thread_potential_energy = Arc::new(RwLock::new(thread_potential_energy));

        (0..atoms.n_atoms).into_par_iter().for_each(|i| {
            let mut tf = thread_forces.write();
            let mut tpe = thread_potential_energy.write();

            for &j in &neighbour_list[i] {
                let mut rij = atoms.positions.column(j) - atoms.positions.column(i);
                atoms.sim_box.apply_boundary_conditions_dis(&mut rij);

                let potential = match self.get_potential_ij(atoms, i, j) {
                    Some(pot) => pot,
                    None => {
                        println!("During force calculation between {} and {} atoms, potential was missing", i + 1, j + 1);
                        continue;
                    }
                };

                let (uij, force_ij) = potential.compute_potential(&rij);
                tpe[i] += uij;
                tf.column_mut(i).add_assign(-force_ij);
            }
        });

        let tf = thread_forces.read();
        atoms.forces.copy_from(&*tf);

        let tpe = thread_potential_energy.read();
        tpe.sum() / 2.0
    }
}

impl_pair_potential_manager!(LJVPBuildListManager);
