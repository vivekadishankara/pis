use na::Matrix3;

use crate::{atoms::new::Atoms, potentials::potential::PotentialManager};

pub enum VelocityDistribution {
    Uniform,
    Gaussian,
}

pub struct StartVelocity {
    pub group: String,
    pub start_velocity: bool,
    pub start_temperature: Option<f64>,
    pub seed: Option<usize>,
    pub dist: Option<VelocityDistribution>,
}

impl Default for StartVelocity {
    fn default() -> Self {
        Self {
            group: String::from("all"),
            start_velocity: true,
            start_temperature: None,
            seed: None,
            dist: None,
        }
    }
}

pub struct NHThermostatChainArgs {
    pub name: String,
    pub group: String,
    pub start_temperature: f64,
    #[allow(dead_code)]
    pub end_temperature: f64,
    pub tau: f64,
}

pub struct MTKBarostatArgs {
    pub name: String,
    pub group: String,
    pub start_pressure: Matrix3<f64>,
    #[allow(dead_code)]
    pub end_pressure: Matrix3<f64>,
    pub tau: f64,
}

pub struct PotentialArgs {
    pub pair_style_args: Vec<String>,
    pub pair_coeff_args: Vec<Vec<String>>,
}

impl Default for PotentialArgs {
    fn default() -> Self {
        let pair_style_args = Vec::new();
        let pair_coeff_args = Vec::new();

        Self {
            pair_style_args,
            pair_coeff_args,
        }
    }
}

pub struct SimulationContext {
    pub atoms: Option<Atoms>,
    pub timestep: f64,
    pub mgr: Option<Box<dyn PotentialManager>>,
    pub starting_velocity: Option<StartVelocity>,
    pub steps: usize,
    pub potential_args: Option<PotentialArgs>,
    pub nh_chain_args: Option<NHThermostatChainArgs>,
    pub mtk_barostat_args: Option<MTKBarostatArgs>,
    // pub dump_style: DumpStyle,
}

impl Default for SimulationContext {
    fn default() -> Self {
        Self {
            atoms: None,
            timestep: 0.01,
            mgr: None,
            starting_velocity: None,
            steps: 100,
            potential_args: None,
            nh_chain_args: None,
            mtk_barostat_args: None,
        }
    }
}

impl SimulationContext {
    pub fn run(&mut self) {
        let mgr = self.mgr.take().expect("Potential Manager not initialized");
        mgr.run(self);
        self.mgr = Some(mgr);
    }
}
