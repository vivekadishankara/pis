use na::Matrix3;

use crate::{
    atoms::new::Atoms,
    constants::KB_KJPERMOLEKELVIN,
    math::symmetrize,
    readers::simulation_context::{MTKBarostatArgs, NHThermostatChainArgs},
    ensemble::nvt::NHThermostatChain,
};

// Martyna, Tobias, Klein (1994) "Constant pressure molecular dynamics algorithms". J. Chem. Phys..
pub struct MTKBarostat {
    #[allow(dead_code)]
    pub name: String,
    #[allow(dead_code)]
    pub group: String,
    pub target_pressure: Matrix3<f64>,
    // barostat velocity
    pub velocity: Matrix3<f64>,
    // barostat mass
    pub w: f64,
    thermostat_chain: NHThermostatChain,
}

impl MTKBarostat {
    pub fn new(
        name: String,
        group: String,
        target_pressure: Matrix3<f64>,
        tau: f64,
        n_atoms: usize,
        thermostat_chain: NHThermostatChain,
    ) -> Self {
        let velocity = Matrix3::zeros();
        let w = ((n_atoms + 1) as f64) * KB_KJPERMOLEKELVIN * thermostat_chain.target_temperature * tau.powi(2);

        Self {
            name,
            group,
            target_pressure,
            velocity,
            w,
            thermostat_chain
        }
    }

    pub fn update_velocity(&mut self, atoms: &Atoms, dt: f64) {
        let instant_pressure = atoms.pressure_tensor();
        let pressure_force =
            (instant_pressure - self.target_pressure) * (atoms.sim_box.volume()) / self.w;
        let mtk_correction = &atoms.kinetic_tensor().diagonal() / atoms.n_atoms as f64;
        let mtk_correction = Matrix3::from_diagonal(&mtk_correction) / self.w;
        let delta_velocity = (pressure_force + mtk_correction) * 0.5 * dt;
        self.velocity += symmetrize(&delta_velocity);
    }

    pub fn scale_h(&self, dt: f64) -> Matrix3<f64> {
        let eta_dot_symmetric = symmetrize(&self.velocity);
        (eta_dot_symmetric * 0.5 * dt).exp()
    }

    pub fn scale_v(&self, dt: f64, particle_n_dof: usize) -> Matrix3<f64> {
        let mut eta_dot_symmetric = symmetrize(&self.velocity);
        let mtk_term2 = (self.velocity.trace() / (particle_n_dof) as f64) * Matrix3::identity();
        eta_dot_symmetric = symmetrize(&(eta_dot_symmetric + mtk_term2));
        (eta_dot_symmetric * -0.5 * dt).exp()
    }

    pub fn kinetic_energy(&self) -> f64 {
        self.w * (self.velocity * self.velocity.transpose()).trace() / 2.0
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
                NHThermostatChain::new_from_args(nh_chain_args, n_atoms).unwrap_or_else(|| {
                    NHThermostatChain::new(
                        "barostat_thermostat".to_string(),
                        "all".to_string(),
                        target_temperature,
                        target_temperature,
                        target_temperature,
                        args.tau / 10.0,
                        3,
                        n_atoms,
                    )
                }),
            )),
            None => None,
        }
    }

    pub fn update_chain(&mut self, dt: f64) {
        // --- First NHC half-step (dt/2) ---
        // Symmetric Trotter:  xi(dt/4) → v(dt/2) → eta(dt/2) → xi(dt/4)
        let n_dof = 6; // 6 because the barostat has 6 degrees of freedom (3 for scaling and 3 for shear)
        let kinetic_energy = self.kinetic_energy();
        self.thermostat_chain.compute_forces(kinetic_energy, n_dof);

        self.thermostat_chain.propagate_xi_backward(0.25 * dt);
        let scale = (-0.5 * dt * self.thermostat_chain.xi[0]).exp();
        self.velocity = &self.velocity * scale;
        self.thermostat_chain.propagate_eta(0.5 * dt);

        let kinetic_energy = self.kinetic_energy();
        self.thermostat_chain.compute_forces(kinetic_energy, n_dof);
        self.thermostat_chain.propagate_xi_forward(0.25 * dt);
    }
}
