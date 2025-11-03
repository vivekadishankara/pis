use na::Matrix3;

use crate::{
    atoms::new::Atoms,
    constants::KB_KJPERMOLEKELVIN,
    math::symmetrize,
    readers::simulation_context::{MTKBarostatArgs, NHThermostatChainArgs},
};

// Martyna, Tobias, Klein (1994) "Constant pressure molecular dynamics algorithms". J. Chem. Phys..
pub struct MTKBarostat {
    #[allow(dead_code)]
    pub name: String,
    #[allow(dead_code)]
    pub group: String,
    pub target_pressure: Matrix3<f64>,
    // barostat momentum
    pub momentum: Matrix3<f64>,
    // barostat mass
    pub w: f64,
}

impl MTKBarostat {
    pub fn new(
        name: String,
        group: String,
        target_pressure: Matrix3<f64>,
        tau: f64,
        n_atoms: usize,
        target_temp: f64,
    ) -> Self {
        let momentum = Matrix3::zeros();
        let w = ((3 * n_atoms) as f64) * KB_KJPERMOLEKELVIN * target_temp * tau.powi(2);

        Self {
            name,
            group,
            target_pressure,
            momentum,
            w,
        }
    }

    pub fn delta_momentum(&self, atoms: &Atoms, dt: f64) -> Matrix3<f64> {
        let instant_pressure = atoms.pressure_tensor();
        let delta_momentum =
            (instant_pressure - self.target_pressure) * (atoms.sim_box.volume() * 0.5 * dt);
        symmetrize(&delta_momentum)
    }

    pub fn scale(&self, dt: f64, velocity_scaling: bool) -> Matrix3<f64> {
        let mut eta_dot = self.momentum / self.w;
        eta_dot = symmetrize(&eta_dot);
        // let scale = mat_exp_taylor(&(-eta_dot * 0.5 * dt));
        let factor = if velocity_scaling { -0.5 } else { 1.0 };
        (eta_dot * factor * dt).exp()
    }

    pub fn kinetic_energy(&self) -> f64 {
        (self.momentum * self.momentum.transpose()).trace() / (2.0 * self.w)
    }

    pub fn potential_energy(&self, h: &Matrix3<f64>) -> f64 {
        (self.target_pressure.transpose() * h).trace()
    }

    pub fn new_from_args(
        mtk_barostat_args: &Option<MTKBarostatArgs>,
        nh_chain_args: &Option<NHThermostatChainArgs>,
        n_atoms: usize,
    ) -> Option<Self> {
        let target_temperature = match nh_chain_args {
            Some(args) => args.start_temperature,
            None => 300.0,
        };
        match mtk_barostat_args {
            Some(args) => Some(Self::new(
                args.name.clone(),
                args.group.clone(),
                args.start_pressure.clone(),
                args.tau,
                n_atoms,
                target_temperature,
            )),
            None => None,
        }
    }
}
