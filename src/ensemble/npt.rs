use na::Matrix3;

use crate::constants::KB_KJPERMOLEKELVIN;

// Martyna, Tobias, Klein (1994) "Constant pressure molecular dynamics algorithms". J. Chem. Phys..
pub struct MTKBarostat {
    pub target_pressure: Matrix3<f64>,
    // barostat coordinate for volume scaling
    pub eta: Matrix3<f64>,
    // barostat momentum
    pub momentum: Matrix3<f64>,
    // barostat mass
    pub w: f64,
}

impl MTKBarostat {
    pub fn new(target_pressure: Matrix3<f64>, tau: f64, n_atoms: usize, target_temp: f64) -> Self {
        let eta = Matrix3::zeros();
        let momentum = Matrix3::zeros();
        let w = ((3 * n_atoms) as f64) * KB_KJPERMOLEKELVIN * target_temp * tau.powi(2);

        Self {
            target_pressure,
            eta,
            momentum,
            w,
        }
    }
}