use na::{DVector, Matrix3xX};

use crate::simulation_box::SimulationBox;

pub struct Atoms {
    pub n_atoms: usize,
    pub type_ids: DVector<usize>,
    pub masses: Vec<f64>,
    pub positions: Matrix3xX<f64>,
    pub velocities: Matrix3xX<f64>,
    pub forces: Matrix3xX<f64>,
    pub sim_box: SimulationBox,
}
