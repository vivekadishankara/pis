use na::Matrix3;

use crate::atoms::new::Atoms;

impl Atoms {
    pub fn scale_box(&mut self, scale: &Matrix3<f64>) {
        let s = &self.sim_box.h_inv * &self.positions;
        self.sim_box.h = scale * self.sim_box.h;
        self.sim_box.h_inv = self
            .sim_box
            .h
            .try_inverse()
            .expect("Box matrix should be invertible");
        self.positions = self.sim_box.h * s;
    }
}
