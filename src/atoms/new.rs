use na::{DVector, Matrix3xX, Vector3};

use crate::simulation_box::SimulationBox;

pub struct Atoms {
    pub n_atoms: usize,
    pub n_types: usize,
    pub positions: Matrix3xX<f64>,
    pub velocities: Matrix3xX<f64>,
    pub forces: Matrix3xX<f64>,
    pub masses: DVector<f64>,
    pub types: Vec<String>,
    pub sim_box: SimulationBox,
}

impl Atoms {
    pub fn default(n_atoms: usize) -> Self {
        Self {
            n_atoms: n_atoms,
            n_types: 1,
            positions: Matrix3xX::zeros(n_atoms),
            velocities: Matrix3xX::zeros(n_atoms),
            forces: Matrix3xX::zeros(n_atoms),
            masses: DVector::from_element(n_atoms, 1.0),
            types: (0..n_atoms).map(|_| String::with_capacity(3)).collect(),
            sim_box: SimulationBox::default(),
        }
    }

    pub fn new(temperature: f64) -> Self {
        let n_atoms: usize = 10;
        let mut this = Self::default(n_atoms);
        this.start_velocities(temperature);
        this
    }
}
