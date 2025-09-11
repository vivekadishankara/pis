use na::{Const, Matrix3, MatrixViewMut, Vector3};

pub struct SimulationBox {
    pub h: Matrix3<f64>,
    h_inv: Matrix3<f64>,
    pub pbc: [bool; 3],
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

    pub fn apply_boundary_conditions_dis(&self, r_ij: &mut Vector3<f64>) {
        let mut s = self.h_inv * *r_ij;
        
        for i in 0..3 {
            if self.pbc[i] {
                s[i] -= s[i].round();
            }
        }

        *r_ij = self.h * s;
    }

    pub fn apply_boundary_conditions_pos<'a>(
        &self, 
        mut r_i: MatrixViewMut<'a, f64, Const<3>, Const<1>>
    ) {
        let mut s = self.h_inv * &r_i;
        
        for i in 0..3 {
            if self.pbc[i] {
                s[i] -= s[i].floor();
            }
        }

        r_i.copy_from(&(self.h * s));
    }

    pub fn from_lammps_data(
        xlo: f64, xhi: f64, 
        ylo: f64, yhi: f64, 
        zlo: f64, zhi: f64, 
        xy: f64, xz: f64, yz: f64
    ) -> Self {
        let ax = xhi - xlo;
        let ay = yhi - ylo;
        let az = zhi - zlo;

        let a = Vector3::new(ax, 0.0, 0.0);
        let b = Vector3::new(xy, ay, 0.0);
        let c = Vector3::new(xz, yz, az);

        let h = Matrix3::from_columns(&[a, b, c]);
        Self::new(h, [true; 3])
    }
}
