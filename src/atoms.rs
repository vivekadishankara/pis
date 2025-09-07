use na::Matrix3xX;

pub struct Atoms {
    n_atoms: usize,
    positions: Matrix3xX<f64>,
    velocities: Matrix3xX<f64>,
    forces: Matrix3xX<f64>,
    masses: Vec<f64>,
    types: Vec<String>,
}

impl Atoms {
    pub fn new_zeroes(n_atoms: usize) -> Self {
        Self {
            n_atoms: n_atoms,
            positions: Matrix3xX::zeros(n_atoms),
            velocities: Matrix3xX::zeros(n_atoms),
            forces: Matrix3xX::zeros(n_atoms),
            masses: vec![1.0; n_atoms],
            types: (0..n_atoms).map(|_| String::with_capacity(3)).collect(),
        }
    }

    pub fn new() -> Self {
        let n_atoms: usize = 10;
        Self::new_zeroes(n_atoms)
    }
}