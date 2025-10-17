use std::{collections::HashSet, ops::AddAssign, sync::Arc};
use na::{DVector, Matrix3xX};
use parking_lot::RwLock;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::{atoms::{neighbour_list::FORWARD_NEIGHBOUR_OFFSETS, new::Atoms}, writers::dump_traj::DumpTraj};
use crate::potentials::potential::PairPotentialManager;


impl Atoms {
    pub fn compute_potential_neighbour_list_parallel(&mut self) -> f64 {
        let neighbour_list = self.build_neighbour_list();

        let thread_forces: Matrix3xX<f64> = Matrix3xX::zeros(self.n_atoms);
        let thread_forces = Arc::new(RwLock::new(thread_forces));

        let thread_potential_energy: DVector<f64> = DVector::zeros(self.n_atoms);
        let thread_potential_energy = Arc::new(RwLock::new(thread_potential_energy));

        (0..self.n_atoms).into_par_iter().for_each(|i| {
            let mut tf = thread_forces.write();
            let mut tpe = thread_potential_energy.write();

            for &j in &neighbour_list[i] {
                let mut rij = self.positions.column(j) - self.positions.column(i);
                self.sim_box.apply_boundary_conditions_dis(&mut rij);

                let potential = match self.get_potential_ij(i, j) {
                    Some(pot) => pot,
                    None => {
                        println!("During force calculation between {} and {} atoms, potential was missing", i + 1, j + 1);
                        continue;
                    }
                };

                let (uij, force_ij) = potential.compute_potetial(&rij);
                tpe[i] += uij;
                tf.column_mut(i).add_assign(-force_ij);
            }
        });

        let tf = thread_forces.read();
        self.forces.copy_from(&*tf);

        let tpe = thread_potential_energy.read();
        tpe.sum() / 2.0
    }

    pub fn build_neighbour_list(&self) -> Vec<Vec<usize>> {
        let max_rcut = self.potential_manager.max_rcut();
        let (nx, ny, nz) = self.divide_into_cells(max_rcut);
        let total_cells = nx * ny * nz;
        let shift = [-1, 0, 1];

        let cells = self.rcut_cells(nx, ny, nz);

        let cell_indices = Self::cell_indices_for_parallel(nx, ny, nz);

        let n_threads = rayon::current_num_threads();
        let neighbour_list_len = (2 * self.n_atoms) / total_cells;
        let thread_neighbour_list: Vec<Vec<usize>> = (0..self.n_atoms)
            .map(|_| Vec::with_capacity(neighbour_list_len))
            .collect();
        let thread_neighbour_list = std::sync::Arc::new(parking_lot::RwLock::new(thread_neighbour_list));

        let thread_neighbour_cells: Vec<HashSet<usize>> = (0..n_threads)
            .map(|_| HashSet::new())
            .collect();
        let thread_neighbour_cells = std::sync::Arc::new(parking_lot::RwLock::new(thread_neighbour_cells));

        cell_indices.par_iter().for_each(|&(cx_i, cy_i, cz_i)| {
            let tid = rayon::current_thread_index().unwrap();

            let mut tnc = thread_neighbour_cells.write();
            let mut tnl = thread_neighbour_list.write();

            let current_cell = Self::cell_index(cx_i, cy_i, cz_i, nx, ny);
            let current_cell_atoms = &cells[current_cell];

            for &dx in &shift {
                for &dy in &shift {
                    for &dz in &shift {
                        let cx_j = (cx_i as isize + dx).rem_euclid(nx as isize) as usize;
                        let cy_j = (cy_i as isize + dy).rem_euclid(ny as isize) as usize;
                        let cz_j = (cz_i as isize + dz).rem_euclid(nz as isize) as usize;

                        let a_neighbour_cell = Self::cell_index(cx_j, cy_j, cz_j, nx, ny);

                        if tnc[tid].contains(&a_neighbour_cell) { continue; }
                        tnc[tid].insert(a_neighbour_cell);
                        let neighbour_cell_atoms = &cells[a_neighbour_cell];

                        for &i in current_cell_atoms.iter() {
                            for &j in neighbour_cell_atoms.iter() {
                                if current_cell == a_neighbour_cell && i==j { continue; }
                                let mut rij = self.positions.column(j) - self.positions.column(i);
                                self.sim_box.apply_boundary_conditions_dis(&mut rij);

                                let potential = match self.get_potential_ij(i, j) {
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

    pub fn compute_potential_verlet_list_parallel(&mut self) -> f64 {
        let max_rcut = self.potential_manager.max_rcut();
        let (nx, ny, nz) = self.divide_into_cells(max_rcut);
        let cells = self.rcut_cells(nx, ny, nz);

        let cell_indices = Self::cell_indices_for_parallel(nx, ny, nz);

        let n_threads = rayon::current_num_threads();

        let thread_forces: Vec<Matrix3xX<f64>> = (0..n_threads)
            .map(|_| Matrix3xX::zeros(self.n_atoms))
            .collect();
        let thread_forces = std::sync::Arc::new(parking_lot::RwLock::new(thread_forces));

        let thread_potential_energy: DVector<f64> = DVector::zeros(n_threads);
        let thread_potential_energy = std::sync::Arc::new(parking_lot::RwLock::new(thread_potential_energy));

        let thread_neighbour_cells: Vec<HashSet<usize>> = (0..n_threads)
            .map(|_| HashSet::new())
            .collect();
        let thread_neighbour_cells = std::sync::Arc::new(parking_lot::RwLock::new(thread_neighbour_cells));

        cell_indices.par_iter().for_each(|&(cx_i, cy_i, cz_i)| {
            let tid = rayon::current_thread_index().unwrap();
            let mut tf = thread_forces.write();
            let mut tpe = thread_potential_energy.write();
            let mut tnc = thread_neighbour_cells.write();

            let current_cell = Self::cell_index(cx_i, cy_i, cz_i, nx, ny);
            let current_cell_atoms = &cells[current_cell];

            for offset in FORWARD_NEIGHBOUR_OFFSETS.iter() {
                let cx_j = (cx_i as isize + offset.dx).rem_euclid(nx as isize) as usize;
                let cy_j = (cy_i as isize + offset.dy).rem_euclid(ny as isize) as usize;
                let cz_j = (cz_i as isize + offset.dz).rem_euclid(nz as isize) as usize;

                let a_neighbour_cell = Self::cell_index(cx_j, cy_j, cz_j, nx, ny);
                if tnc[tid].contains(&a_neighbour_cell) { continue; }
                tnc[tid].insert(a_neighbour_cell);
                let neighbour_cell_atoms = &cells[a_neighbour_cell];

                for &i in current_cell_atoms.iter() {
                    for &j in neighbour_cell_atoms.iter() {
                        if current_cell == a_neighbour_cell && i<=j { continue; }
                        let mut rij = self.positions.column(j) - self.positions.column(i);
                        self.sim_box.apply_boundary_conditions_dis(&mut rij);

                        let potential = match self.get_potential_ij(i, j) {
                            Some(pot) => pot,
                            None => {
                                println!("During force calculation between {} and {} atoms, potential was missing", i + 1, j + 1);
                                continue;
                            }
                        };
                        
                        if rij.norm() > potential.get_rcut() { continue; }

                        let (uij, force_ij) = potential.compute_potetial(&rij);

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
            self.forces += &tf[i];
        }
        
        let tpe = thread_potential_energy.read();
        tpe.sum()

    }
}
