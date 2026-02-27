//! This module defines the struct Atoms
use na::{DVector, Matrix3xX};

use crate::simulation_box::SimulationBox;

/// This struct holds the atoms in the simulation.
/// Atom masses, positions, velocities and other defining characteristics are the member
/// variables of this struct.
pub struct Atoms {
    pub n_atoms: usize,
    pub type_ids: DVector<usize>,
    pub masses: Vec<f64>,
    pub positions: Matrix3xX<f64>,
    pub velocities: Matrix3xX<f64>,
    pub forces: Matrix3xX<f64>,
    pub sim_box: SimulationBox,
}
