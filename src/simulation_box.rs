use na::{Matrix3, Vector3};

pub struct SimulationBox {
    h: Matrix3<f64>,
    h_inv: Matrix3<f64>,
    pbc: [bool; 3],
}

impl SimulationBox {
    pub fn new(h: Matrix3<f64>, pbc: [bool; 3]) -> Self {
        let h_inv = h.try_inverse().expect("Box matrix should be invertible");
        Self {
            h,
            h_inv,
            pbc,
        }
    }

    pub fn default() -> Self {
        Self {
            h: Matrix3::identity(),
            h_inv: Matrix3::identity(),
            pbc: [true; 3],
        }
    }

    pub fn apply_boundary_conditions(&self, rij: &Vector3<f64>) -> Vector3<f64> {
        let mut s = self.h_inv * rij;

        for i in 0..3 {
            if self.pbc[i] {
                s[i] -= s[i].round();
            }
        }

        self.h * s
    }
}