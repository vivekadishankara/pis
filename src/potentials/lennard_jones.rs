use na::Vector3;

pub struct Potentials;

impl Potentials {
    pub fn lennard_jones(
        rij: &Vector3<f64>, 
        sigma: f64, 
        epsilon: f64,
        rcut: f64,
        shift: bool,
    ) -> (f64, Vector3<f64>){
        let rij2 = rij.norm_squared();
        let inv_rij2 = 1.0 / rij2;
        let vanderwaals_attraction = (sigma.powi(2) * inv_rij2).powi(3);
        let lj_repulsion = vanderwaals_attraction.powi(2);

        let mut potential_energy = 4.0 * epsilon * (lj_repulsion - vanderwaals_attraction);

        let force = 24.0 * epsilon * (2.0 * lj_repulsion - vanderwaals_attraction) * inv_rij2 * rij;

        if shift {
            let cutoff_inv2 = (sigma / rcut).powi(2);
            let cutoff_attraction = cutoff_inv2.powi(3);
            let cutoff_repulsion = cutoff_attraction.powi(2);

            let u_cutoff = 4.0 * epsilon * (cutoff_repulsion - cutoff_attraction);

            potential_energy -= u_cutoff;
        }

        (potential_energy, force)
    }
}
