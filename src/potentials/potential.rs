use core::panic;
use na::{Matrix3xX, Vector3};
use std::collections::HashMap;

use crate::atoms::new::Atoms;
use crate::ensemble::npt::MTKBarostat;
use crate::ensemble::nvt::NHThermostatChain;
use crate::math::{symmetrize, mat_exp_taylor};
use crate::writers::dump_traj::DumpTraj;

pub trait PotentialManager: Send + Sync {
    fn compute_potential(&self, atoms: &mut Atoms) -> f64;

    fn verlet_step_nve(&self, atoms: &mut Atoms, dt: f64) -> f64 {
        let a_t = atoms.current_acceleration();

        atoms.positions += &atoms.velocities * dt + &a_t * 0.5 * dt.powi(2);

        for r_i in atoms.positions.column_iter_mut() {
            atoms.sim_box.apply_boundary_conditions_pos(r_i);
        }

        atoms.forces = Matrix3xX::zeros(atoms.n_atoms);

        let potential_energy = self.compute_potential(atoms);

        let a_tdt = atoms.current_acceleration();

        atoms.velocities += (a_t + a_tdt) * 0.5 * dt;

        potential_energy
    }

    fn verlet_step_nvt_nhc(
        &self,
        atoms: &mut Atoms,
        dt: f64,
        noose_hoover_chain: &mut NHThermostatChain,
    ) -> f64 {
        let mut kinetic_energy = atoms.kinetic_energy();

        noose_hoover_chain.compute_forces(kinetic_energy, atoms.n_atoms);
        noose_hoover_chain.propagate_half_step(dt);
        let scale = (-0.5 * dt * noose_hoover_chain.xi[0]).exp();
        atoms.velocities = &atoms.velocities * scale;

        let potential_energy = self.verlet_step_nve(atoms, dt);

        atoms.velocities = &atoms.velocities * scale;

        kinetic_energy = atoms.kinetic_energy();

        noose_hoover_chain.compute_forces(kinetic_energy, atoms.n_atoms);
        noose_hoover_chain.propagate_half_step(dt);

        potential_energy
    }

    fn verlet_step_npt_mtk(&self, atoms: &mut Atoms, dt: f64, mtk_barostat: &mut MTKBarostat, noose_hoover_chain: &mut NHThermostatChain) {
        let mut kinetic_energy = atoms.kinetic_energy();
        let mut instant_pressure = atoms.pressure_tensor();

        let delta_momentum = atoms.sim_box.volume() * (instant_pressure - mtk_barostat.target_pressure) * 0.5 * dt;
        mtk_barostat.momentum += symmetrize(&delta_momentum);

        let mut eta_dot = mtk_barostat.momentum / mtk_barostat.w;
        eta_dot = symmetrize(&eta_dot);
        let mut scale = mat_exp_taylor(&(-eta_dot * 0.5 * dt));

        atoms.velocities = &scale * &atoms.velocities;

        self.verlet_step_nvt_nhc(atoms, dt, noose_hoover_chain);

        mtk_barostat.eta += eta_dot * dt;
        scale = (-eta_dot * 0.5 * dt / 3.0).exp();

        atoms.velocities = &atoms.velocities * scale;

        kinetic_energy = atoms.kinetic_energy();
        instant_pressure = atoms.pressure_tensor();

        mtk_barostat.momentum += atoms.sim_box.volume() * (instant_pressure - mtk_barostat.target_pressure) * 0.5 * dt;
    }

    fn run(&self, atoms: &mut Atoms, dt: f64, time_steps: usize, dump_path: &str) {
        let mut dumper = DumpTraj::new(dump_path).expect("Failed to create dump file");
        dumper.write_step(&atoms, 0).expect("Failed to write step");
        let first_potential = self.compute_potential(atoms);
        println!("{} {}", 0, first_potential);

        let ensemble = "nvt";

        let mut noose_hoover_chain = NHThermostatChain::new(5.0, 100.0, 3);

        for i in 0..time_steps {
            let step_potential = match ensemble {
                "nve" => self.verlet_step_nve(atoms, dt),
                "nvt" => self.verlet_step_nvt_nhc(atoms, dt, &mut noose_hoover_chain),
                _ => panic!("Ensemble unknown"),
            };
            dumper
                .write_step(&atoms, i + 1)
                .expect("Failed to write step");

            let kinetic_energy = atoms.kinetic_energy();
            let basic_hamiltonian = step_potential + kinetic_energy;
            let hamiltonian = match ensemble {
                "nve" => basic_hamiltonian,
                "nvt" => {
                    basic_hamiltonian
                        + noose_hoover_chain.kinetic_energy()
                        + noose_hoover_chain.potential_energy(atoms.n_atoms)
                }
                _ => panic!("Ensemble unknown"),
            };
            let temperature = atoms.temerature(kinetic_energy);
            let pressure = atoms.pressure(kinetic_energy);

            println!(
                "{} {:.3} {:.3} {:.3} {:.3} {:.3}",
                i + 1,
                step_potential,
                kinetic_energy,
                hamiltonian,
                temperature,
                pressure
            );
        }
    }
}

pub trait PairPotential: Send + Sync {
    fn compute_potential(&self, rij: &Vector3<f64>) -> (f64, Vector3<f64>);
    fn get_rcut(&self) -> f64;
}

type AtomPair = (usize, usize);
pub type Table = HashMap<AtomPair, Box<dyn PairPotential>>;

pub trait PairPotentialManager: Sized {
    fn with_table(table: Table) -> Self;
    fn table(&self) -> &Table;
    fn table_mut(&mut self) -> &mut Table;

    fn new() -> Self {
        Self::with_table(Table::new())
    }

    fn insert<P>(&mut self, key: AtomPair, potential: P)
    where
        P: PairPotential + 'static,
    {
        self.table_mut().insert(key, Box::new(potential));
    }

    fn get(&self, key: &AtomPair) -> Option<&dyn PairPotential> {
        self.table().get(key).map(|b| b.as_ref())
    }

    #[allow(dead_code)]
    fn max_rcut(&self) -> f64 {
        let mut max_rcut = 0.0;
        for potential in self.table().values() {
            let rcut = (*potential).get_rcut();
            if max_rcut < rcut {
                max_rcut = rcut;
            }
        }
        max_rcut
    }

    fn get_potential_ij(&self, atoms: &Atoms, i: usize, j: usize) -> Option<&dyn PairPotential> {
        let type_i = atoms.type_ids[i];
        let type_j = atoms.type_ids[j];

        let type_tuple = if type_i < type_j {
            (type_i, type_j)
        } else {
            (type_j, type_i)
        };

        self.get(&type_tuple)
    }
}
