use na::Matrix3;

use crate::{atoms::new::Atoms, potentials::potential::PotentialManager};

/// An enum for the style of velocity distribution
pub enum VelocityDistribution {
    Uniform,
    Gaussian,
}

/// Context for the "velocity" command
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

/// Context for the thermostat for the fix command for nvt and npt style
pub struct NHThermostatChainArgs {
    pub name: String,
    pub group: String,
    pub start_temperature: f64,
    #[allow(dead_code)]
    pub end_temperature: f64,
    pub tau: f64,
}

/// Context for the barostat for the fix command for npt style
pub struct MTKBarostatArgs {
    pub name: String,
    pub group: String,
    pub start_pressure: Matrix3<f64>,
    #[allow(dead_code)]
    pub end_pressure: Matrix3<f64>,
    pub tau: f64,
}

/// Context for the Potential arguments combining the arguments for "pair_style" and "pair_coeff" command
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

/// Context for the "dump" command
pub struct DumpArgs {
    pub name: String,
    pub group: String,
    pub style: String,
    pub dump_step: usize,
    pub file_name: String,
}

impl Default for DumpArgs {
    fn default() -> Self {
        let name = String::from("default_dump");
        let group = String::from("all");
        let style = String::from("atoms");
        let dump_step = 1;
        let file_name = String::from("dump.lammpstrj");

        Self {
            name,
            group,
            style,
            dump_step,
            file_name,
        }
    }
}

/// Context for the entire molecular dynamic run
pub struct SimulationContext {
    pub atoms: Option<Atoms>,
    pub timestep: f64,
    pub mgr: Option<Box<dyn PotentialManager>>,
    pub starting_velocity: Option<StartVelocity>,
    pub steps: usize,
    pub potential_args: Option<PotentialArgs>,
    pub nh_chain_args: Option<NHThermostatChainArgs>,
    pub mtk_barostat_args: Option<MTKBarostatArgs>,
    pub dump_args: DumpArgs,
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
            dump_args: DumpArgs::default(),
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
