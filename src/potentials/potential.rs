use na::{Matrix3xX, Vector3};
use std::collections::HashMap;

use crate::atoms::new::Atoms;
use crate::ensemble::npt::MTKBarostat;
use crate::ensemble::nvt::NHThermostatChain;
use crate::potentials::kind::PairPotentialKind;

/// General interface for anything that can compute forces and energies.
/// For pair potentials specifically, also implement [`PairPotentialManager`]
/// which provides the pair table infrastructure.
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
        // --- First NHC half-step (dt/2) ---
        // Symmetric Trotter:  xi(dt/4) → v(dt/2) → eta(dt/2) → xi(dt/4)
        noose_hoover_chain.half_step(atoms, dt);

        // --- NVE step (full dt) ---
        let potential_energy = self.verlet_step_nve(atoms, dt);

        // --- Second NHC half-step (dt/2) ---
        noose_hoover_chain.half_step(atoms, dt);

        potential_energy
    }

    fn verlet_step_npt_mtk(
        &self,
        atoms: &mut Atoms,
        dt: f64,
        mtk_barostat: &mut MTKBarostat,
        noose_hoover_chain: &mut NHThermostatChain,
    ) -> f64 {
        // 1. Advance Thermostats (Both particle and barostat NH chains) by dt/2
        noose_hoover_chain.half_step(atoms, dt);
        mtk_barostat.update_chain(dt);

        // 2. Advance Barostat velocity (driving pressure) by dt/2
        mtk_barostat.update_velocity(atoms, dt);

        // 3. Scale particle velocities by the MTK factor for dt/2
        let scale = mtk_barostat.scale_v(dt, atoms.degress_of_freedom());
        atoms.velocities = &scale * &atoms.velocities;

        // 4. Standard Kick: Advance velocities using current forces for dt/2
        let a_t = atoms.current_acceleration();
        atoms.velocities += &a_t * 0.5 * dt;

        // 5. Symmetric Drift: Advance Box and Positions simultaneously
        let scale_h = mtk_barostat.scale_h(dt);
        atoms.scale_box(&scale_h);

        // Main drift step
        atoms.positions += &atoms.velocities * dt;

        for r_i in atoms.positions.column_iter_mut() {
            atoms.sim_box.apply_boundary_conditions_pos(r_i);
        }
        
        // Advance the box matrix and positions by the remaining half-step factor
        let scale_h = mtk_barostat.scale_h(dt);
        atoms.scale_box(&scale_h);

        // 6. Force Evaluation at new positions
        atoms.forces = Matrix3xX::zeros(atoms.n_atoms);
        let potential_energy = self.compute_potential(atoms);
        let a_tdt = atoms.current_acceleration();

        // 7. Standard Kick: Advance velocities using NEW forces for dt/2
        atoms.velocities += &a_tdt * 0.5 * dt;

        // 8. Scale particle velocities by the MTK factor for remaining dt/2
        let scale = mtk_barostat.scale_v(dt, atoms.degress_of_freedom());
        atoms.velocities = &scale * &atoms.velocities;

        // 9. Advance Barostat velocity by remaining dt/2
        mtk_barostat.update_velocity(atoms, dt);

        // 10. Clean up Thermostats for the final dt/2
        mtk_barostat.update_chain(dt);
        noose_hoover_chain.half_step(atoms, dt);

        potential_energy
    }
}

pub trait PairPotential: Send + Sync {
    fn compute_potential(&self, rij: &Vector3<f64>) -> (f64, Vector3<f64>);
    fn get_rcut(&self) -> f64;
}

type AtomPair = (usize, usize);
pub type Table = HashMap<AtomPair, PairPotentialKind>;

/// Additional interface for potentials that decompose into pairwise interactions.
/// Implementors should also implement [`PotentialManager`].
pub trait PairPotentialManager: Sized {
    fn with_table(table: Table) -> Self;
    fn table(&self) -> &Table;
    fn table_mut(&mut self) -> &mut Table;

    fn new() -> Self {
        Self::with_table(Table::new())
    }

    fn is_empty(&self) -> bool {
        self.table().is_empty()
    }

    fn insert(&mut self, key: AtomPair, potential: PairPotentialKind) {
        self.table_mut().insert(key, potential);
    }

    fn get(&self, key: &AtomPair) -> Option<&PairPotentialKind> {
        self.table().get(key)
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

    fn get_potential_ij(&self, atoms: &Atoms, i: usize, j: usize) -> Option<&PairPotentialKind> {
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

macro_rules! impl_pair_potential_manager {
    ($name:ident) => {
        impl PairPotentialManager for $name {
            fn with_table(table: Table) -> Self {
                Self { table }
            }
            fn table(&self) -> &Table {
                &self.table
            }
            fn table_mut(&mut self) -> &mut Table {
                &mut self.table
            }
        }
    };
}
