use na::{DVector, Matrix3xX};

pub struct Atoms {
    pub n_atoms: usize,
    pub positions: Matrix3xX<f64>,
    pub velocities: Matrix3xX<f64>,
    pub forces: Matrix3xX<f64>,
    pub masses: DVector<f64>,
    pub types: Vec<String>,
}

impl Atoms {
    pub fn new_zeroes(n_atoms: usize) -> Self {
        Self {
            n_atoms: n_atoms,
            positions: Matrix3xX::zeros(n_atoms),
            velocities: Matrix3xX::zeros(n_atoms),
            forces: Matrix3xX::zeros(n_atoms),
            masses: DVector::from_element(n_atoms, 1.0),
            types: (0..n_atoms).map(|_| String::with_capacity(3)).collect(),
        }
    }

    pub fn new(temperature: f64) -> Self {
        let n_atoms: usize = 10;
        let mut this = Self::new_zeroes(n_atoms);
        this.start_velocities(temperature);
        this
    }
}
