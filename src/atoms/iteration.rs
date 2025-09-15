use std::{collections::HashSet};
use na::{DVector, Matrix3xX};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::{atoms::{neighbour_list::FORWARD_NEIGHBOUR_OFFSETS, new::Atoms}, writers::dump_traj::DumpTraj};


impl Atoms {
    pub fn compute_potential_and_forces_n_list_parallel(&mut self) -> f64 {
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
            let tnc = thread_neighbour_cells.write();

            let current_cell = Self::cell_index(cx_i, cy_i, cz_i, nx, ny);
            let current_cell_atoms = &cells[current_cell];

            for offset in FORWARD_NEIGHBOUR_OFFSETS.iter() {
                let cx_j = (cx_i as isize + offset.dx).rem_euclid(nx as isize) as usize;
                let cy_j = (cy_i as isize + offset.dy).rem_euclid(ny as isize) as usize;
                let cz_j = (cz_i as isize + offset.dz).rem_euclid(nz as isize) as usize;

                let a_neighbour_cell = Self::cell_index(cx_j, cy_j, cz_j, nx, ny);
                if tnc[tid].contains(&a_neighbour_cell) { continue; }
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
        });

        let tf = thread_forces.read();

        for i in 0..n_threads {
            self.forces += &tf[i];
        }
        
        let tpe = thread_potential_energy.read();
        tpe.sum()

    }

    pub fn compute_potential_and_forces_n_list(&mut self) -> f64 {
        let max_rcut = self.potential_manager.max_rcut();
        let (nx, ny, nz) = self.divide_into_cells(max_rcut);
        let cells = self.rcut_cells(nx, ny, nz);

        let mut neighbour_cells: HashSet<usize> = HashSet::new();
        let mut potential_energy: f64 = 0.0;

        for cx_i in 0..nx {
            for cy_i in 0..ny {
                for cz_i in 0..nz {
                    let current_cell = Self::cell_index(cx_i, cy_i, cz_i, nx, ny);
                    let current_cell_atoms = &cells[current_cell];

                    for offset in FORWARD_NEIGHBOUR_OFFSETS.iter() {
                        let cx_j = (cx_i as isize + offset.dx).rem_euclid(nx as isize) as usize;
                        let cy_j = (cy_i as isize + offset.dy).rem_euclid(ny as isize) as usize;
                        let cz_j = (cz_i as isize + offset.dz).rem_euclid(nz as isize) as usize;

                        let a_neighbour_cell = Self::cell_index(cx_j, cy_j, cz_j, nx, ny);
                        if neighbour_cells.contains(&a_neighbour_cell) { continue; }
                        neighbour_cells.insert(a_neighbour_cell);

                        let neighbour_cell_atoms = &cells[a_neighbour_cell];

                        for &i in current_cell_atoms.iter() {
                            for &j in neighbour_cell_atoms.iter() {
                                if current_cell == a_neighbour_cell && i<=j { continue; }
                                let mut rij = self.positions.column(j) - self.positions.column(i);
                                self.sim_box.apply_boundary_conditions_dis(&mut rij);

                                let potential = match self.get_potential_ij(i,j) {
                                    Some(pot) => pot,
                                    None => {
                                        println!("During force calculation between {} and {} atoms, potential was missing", i + 1, j + 1);
                                        continue;
                                    }                                            
                                };

                                if rij.norm() > potential.get_rcut() { continue; }
                                let (uij, force_ij) = potential.compute_potetial(&rij);
                                potential_energy += uij;
                                {
                                    let mut fi = self.forces.column_mut(i);
                                    fi -= force_ij;

                                    let mut fj = self.forces.column_mut(j);
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

    #[allow(dead_code)]
    pub fn compute_potential_and_forces_n_list1 (&mut self) -> f64 {
        let max_rcut = self.potential_manager.max_rcut();

        let (nx, ny, nz) = self.divide_into_cells(max_rcut);
        let shift = [-1, 0, 1];
        let cells = self.rcut_cells(nx, ny, nz);

        let mut neighbour_cells: HashSet<usize> = HashSet::new();
        let mut potential_energy: f64 = 0.0;

        for cx_i in 0..nx {
            for cy_i in 0..ny {
                for cz_i in 0..nz {
                    let current_cell = Self::cell_index(cx_i, cy_i, cz_i, nx, ny);
                    let current_cell_atoms = &cells[current_cell];

                    for dx in &shift {
                        for dy in &shift {
                            for dz in &shift {
                                let cx_j = (cx_i as isize + dx).rem_euclid(nx as isize) as usize;
                                let cy_j = (cy_i as isize + dy).rem_euclid(ny as isize) as usize;
                                let cz_j = (cz_i as isize + dz).rem_euclid(nz as isize) as usize;

                                let a_neighbour_cell = Self::cell_index(cx_j, cy_j, cz_j, nx, ny);
                                if neighbour_cells.contains(&a_neighbour_cell) { continue; }
                                neighbour_cells.insert(a_neighbour_cell);

                                let neighbour_cell_atoms = &cells[a_neighbour_cell];

                                for &i in current_cell_atoms.iter() {
                                    for &j in neighbour_cell_atoms.iter() {
                                        if i <= j { continue; }
                                        let mut rij = self.positions.column(j) - self.positions.column(i);
                                        self.sim_box.apply_boundary_conditions_dis(&mut rij);

                                        let potential = match self.get_potential_ij(i,j) {
                                            Some(pot) => pot,
                                            None => {
                                                println!("During force calculation between {} and {} atoms, potential was missing", i + 1, j + 1);
                                                continue;
                                            }                                            
                                        };

                                        if rij.norm() > potential.get_rcut() { continue; }
                                        let (uij, force_ij) = potential.compute_potetial(&rij);
                                        potential_energy += uij;
                                        {
                                            let mut fi = self.forces.column_mut(i);
                                            fi -= force_ij;

                                            let mut fj = self.forces.column_mut(j);
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

    #[allow(dead_code)]
    pub fn compute_potential_and_forces(&mut self) -> f64 {
        let mut potential_energy: f64 = 0.0;
        for i in 0..self.n_atoms {
            for j in (i + 1)..self.n_atoms {
                let mut rij = self.positions.column(j) - self.positions.column(i);

                self.sim_box.apply_boundary_conditions_dis(&mut rij);
                let potential = match self.get_potential_ij(i, j) {
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
                    let mut fi = self.forces.column_mut(i);
                    fi -= force_ij;

                    let mut fj = self.forces.column_mut(j);
                    fj += force_ij;
                }
            }
        }

        potential_energy
    }

    pub fn verlet_step(&mut self, dt: f64) -> f64 {

        let a_t = self.current_acceleration();

        self.positions += &self.velocities * dt + &a_t * 0.5 * dt.powi(2);

        for r_i in self.positions.column_iter_mut() {
            self.sim_box.apply_boundary_conditions_pos(r_i);
        }
        
        self.forces = Matrix3xX::zeros(self.n_atoms);

        let potential_energy = self.compute_potential_and_forces_n_list_parallel();

        let a_tdt = self.current_acceleration();

        self.velocities += (a_t + a_tdt) * 0.5 * dt;

        potential_energy
    }

    pub fn run(
        &mut self, 
        dt: f64, 
        time_steps: usize,
        dump_path: &str,
    ) -> DVector<f64> {
        let mut potential_energies: DVector<f64> = DVector::zeros(time_steps + 1);
        let mut dumper = DumpTraj::new(dump_path).expect("Failed to create dump file");
        dumper.write_step(self, 0).expect("Failed to write step");
        let first_potential = self.compute_potential_and_forces_n_list_parallel();
        potential_energies[0] = first_potential;
        println!("{} {}", 0, first_potential);
        for i in 0..time_steps {
            let step_potential = self.verlet_step(dt);
            potential_energies[i + 1] = step_potential;
            dumper.write_step(self, i + 1).expect("Failed to write step");
            println!("{} {}", i + 1, step_potential);
        }
        potential_energies
    }
}
