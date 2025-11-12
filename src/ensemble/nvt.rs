use crate::{constants::KB_KJPERMOLEKELVIN, readers::simulation_context::NHThermostatChainArgs};

// implementation according to https://www2.stat.duke.edu/~scs/Projects/REMD/NoseHooverChains1992.pdf
pub struct NHThermostatChain {
    #[allow(dead_code)]
    pub name: String,
    #[allow(dead_code)]
    pub group: String,
    pub chain_size: usize,
    pub start_temperature: f64,
    pub end_temperature: f64,
    pub target_temperature: f64,
    pub xi: Vec<f64>,  // thermostat velocities
    pub eta: Vec<f64>, // thermostat coordinates
    pub g: Vec<f64>,   // thermostat forces
    pub q: Vec<f64>,   // thermostat masses
}

impl NHThermostatChain {
    // tau is the reaxation time for the thermostat chain
    pub fn new(
        name: String,
        group: String,
        start_temperature: f64,
        end_temperature: f64,
        target_temperature: f64,
        tau: f64,
        chain_size: usize,
    ) -> Self {
        // initialise thermostat velocities, positions forces and masses to 0.0
        let xi = vec![0.0; chain_size];
        let eta = vec![0.0; chain_size];
        let g = vec![0.0; chain_size];
        let mut q = vec![0.0; chain_size];

        // mass of the first thermostat
        let q_value = KB_KJPERMOLEKELVIN * target_temperature * tau.powi(2);

        // damp the higher thermostats by a factor of 10
        for i in 0..chain_size {
            q[i] = q_value / (10.0_f64).powi(i as i32);
        }

        Self {
            name,
            group,
            chain_size,
            start_temperature,
            end_temperature,
            target_temperature,
            xi,
            eta,
            g,
            q,
        }
    }

    // Compute generalized thermostat forces
    pub fn compute_forces(&mut self, kinetic_energy: f64, n_atoms: usize) {
        // G1​= 2K − Ndof ​kB ​T
        self.g[0] = 2.0 * kinetic_energy
            - ((n_atoms * 3) as f64) * KB_KJPERMOLEKELVIN * self.target_temperature;

        for j in 1..self.chain_size {
            // Gj ​= Q(j−1) ​ξ(j−1)^2​ − kB​ T for j≥2
            self.g[j] = self.q[j - 1] * self.xi[j - 1].powi(2)
                - KB_KJPERMOLEKELVIN * self.target_temperature;
        }
    }

    // propagation of the thermostats to half a timestep
    pub fn propagate_half_step(&mut self, timestep: f64) {
        let j = self.chain_size - 1;
        self.xi[j] = 0.5 * timestep * self.g[j] / self.q[j];

        for l in (0..j).rev() {
            self.xi[l] = (0.5 * timestep * self.g[l] / self.q[l])
                * (-0.25 * timestep * self.xi[l + 1]).exp();
        }
    }

    pub fn kinetic_energy(&self) -> f64 {
        let mut thermostat_ke = 0.0;
        for i in 0..self.chain_size {
            thermostat_ke += 0.5 * self.q[i] * self.xi[i].powi(2);
        }
        thermostat_ke
    }

    pub fn potential_energy(&self, n_atoms: usize) -> f64 {
        let mut thermostat_pe =
            (n_atoms * 3) as f64 * KB_KJPERMOLEKELVIN * self.target_temperature * self.eta[0];
        for i in 1..self.chain_size {
            thermostat_pe += KB_KJPERMOLEKELVIN * self.target_temperature * self.eta[i];
        }
        thermostat_pe
    }

    pub fn new_from_args(nh_chain_args: &Option<NHThermostatChainArgs>) -> Option<Self> {
        match nh_chain_args {
            Some(args) => Some(Self::new(
                args.name.to_string(),
                args.group.to_string(),
                args.start_temperature,
                args.end_temperature,
                args.start_temperature,
                args.tau,
                3,
            )),
            None => None,
        }
    }

    pub fn calculate_target_temperature(
        &mut self,
        current_timestep: usize,
        total_timesteps: usize,
    ) {
        // T_target = T_start + ((T_end - T_start) / total_timesteps) * current_timestep
        self.target_temperature = self.start_temperature
            + ((self.end_temperature - self.start_temperature) / (total_timesteps) as f64)
                * current_timestep as f64;
    }
}
