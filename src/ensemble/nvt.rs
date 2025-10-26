use crate::constants::KB_KJPERMOLEKELVIN;

// implementation according to https://www2.stat.duke.edu/~scs/Projects/REMD/NoseHooverChains1992.pdf
pub struct NHThermostatChain {
    pub chain_size: usize,
    pub target_temp: f64,
    pub xi: Vec<f64>,  // thermostat velocities
    pub eta: Vec<f64>, // thermostat coordinates
    pub g: Vec<f64>,   // thermostat forces
    pub q: Vec<f64>,   // thermostat masses
}

impl NHThermostatChain {
    // tau is the reaxation time for the thermostat chain
    pub fn new(temperature: f64, tau: f64, chain_size: usize) -> Self {
        // initialise thermostat velocities, positions forces and masses to 0.0
        let xi = vec![0.0; chain_size];
        let eta = vec![0.0; chain_size];
        let g = vec![0.0; chain_size];
        let mut q = vec![0.0; chain_size];

        // mass of the first thermostat
        let q_value = KB_KJPERMOLEKELVIN * temperature * tau.powi(2);

        // damp the higher thermostats by a factor of 10
        for i in 0..chain_size {
            q[i] = q_value / (10.0_f64).powi(i as i32);
        }

        Self {
            chain_size,
            target_temp: temperature,
            xi,
            eta,
            g,
            q,
        }
    }

    // Compute generalized thermostat forces
    pub fn compute_forces(&mut self, kinetic_energy: f64, n_atoms: usize) {
        // G1​= 2K − Ndof ​kB ​T
        self.g[0] =
            2.0 * kinetic_energy - ((n_atoms * 3) as f64) * KB_KJPERMOLEKELVIN * self.target_temp;

        for j in 1..self.chain_size {
            // Gj ​= Q(j−1) ​ξ(j−1)^2​ − kB​ T for j≥2
            self.g[j] =
                self.q[j - 1] * self.xi[j - 1].powi(2) - KB_KJPERMOLEKELVIN * self.target_temp;
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
            (n_atoms * 3) as f64 * KB_KJPERMOLEKELVIN * self.target_temp * self.eta[0];
        for i in 1..self.chain_size {
            thermostat_pe += KB_KJPERMOLEKELVIN * self.target_temp * self.eta[i];
        }
        thermostat_pe
    }
}
