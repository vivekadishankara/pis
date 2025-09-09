use crate::atoms::new::Atoms;

impl Atoms {
    pub fn iterate(&mut self, sigma: f64, epsilon: f64) -> f64 {
        for i in 0..self.n_atoms {
            for j in (i + 1)..self.n_atoms {
                let rij = self.positions.column(i) - self.positions.column(j);

                let pbc_rij = self.sim_box.apply_boundary_conditions(&rij);
            }
        }

        0.0
    }
}