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
    pub dist: Option<VelocityDistribution>
}

impl Default for StartVelocity {
    fn default() -> Self {
        Self { group: String::from("all"), start_velocity: true, start_temperature: None, seed: None, dist: None }
    }
}

pub struct SimulationContext {
    pub atoms: Option<Atoms>,
    pub timestep: f64,
    pub mgr: Option<Box<dyn PotentialManager>>,
    pub starting_velocity: Option<StartVelocity>,
    pub steps: usize,
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
        }
    }
}