use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
};

use crate::{
    potentials::{
        lennard_jones::{LJVOffsetManager, LennardJones},
        potential::PairPotentialManager,
    },
    readers::{
        input_file::commands::{
            Command, Dump, Fix, PairCoeff, PairStyle, ReadData, RunSteps, TimeStep, Velocity,
        },
        simulation_context::SimulationContext,
    },
};


/// [`System`] is the basic API for running the molecular dynamics solution. 
/// 
/// It's [`System::new`] function takes in the path to the input file and initializes the system.
/// We then need to contextualize the system according to the arguments given in the input file. This is the process where the atoms are initialized, the thermostat and barostats are set and so on,
/// The contextualized system then needs to run which is the molecular dynamics run for the sytem created from before.
/// 
/// # Examples
/// 
/// A typical main then looks like this:
/// 
/// ```
/// use crate::system::System;
/// 
/// fn main() {
///     System::new(path_to_file).read().contextualize().run();
/// }
/// ```
pub struct System {
    /// the path to the input file which contains the arguments to initialize and run the molecular system
    infile: String,
    command_hash: HashMap<String, Box<dyn Command>>,
    ctx: SimulationContext,
}

impl System {
    /// The constructor for the System which takes in the path to the input file
    pub fn new(infile: String) -> Self {
        let command_hash = HashMap::new();
        let ctx = SimulationContext::default();
        Self {
            infile,
            command_hash,
            ctx,
        }
    }

    /// This function is where we add the string corresponding to the command in the input file. 
    /// The parsing of arguments is done in a struct which implements the [`Command`] trait. 
    fn build_command_hash(&mut self) {
        self.command_hash
            .insert(String::from("timestep"), Box::new(TimeStep));
        self.command_hash
            .insert(String::from("run"), Box::new(RunSteps));
        self.command_hash
            .insert(String::from("velocity"), Box::new(Velocity));
        self.command_hash
            .insert(String::from("read_data"), Box::new(ReadData));
        self.command_hash
            .insert(String::from("pair_style"), Box::new(PairStyle));
        self.command_hash
            .insert(String::from("pair_coeff"), Box::new(PairCoeff));
        self.command_hash.insert(String::from("fix"), Box::new(Fix));
        self.command_hash
            .insert(String::from("dump"), Box::new(Dump));
    }

    /// Reads the input file and collects all the arguments provided by the input file.
    pub fn read(&mut self) -> &mut Self {
        let file = File::open(&self.infile).expect("Failed to open {infile}");
        let reader = BufReader::new(file);

        self.build_command_hash();

        for line in reader.lines() {
            let line = line.unwrap();
            let line = line.trim();

            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            // Considering only the part of the line that is before the commented part.
            let uncommented = line.splitn(2, '#').next().unwrap().trim();

            let line_split: Vec<&str> = uncommented.split_whitespace().collect();

            let command = line_split[0];
            let args = &line_split[1..];

            if let Some(handler) = self.command_hash.get(command) {
                let _ = handler.run(args, &mut self.ctx);
            } else {
                println!("Command Unkown");
            }
        }
        self
    }

    /// Contextualises the system with respect to the arguments read in the [`System::read`].
    /// It actually creates the molecur system to be run.
    // TODO: Maybe the function can be diluted out since it only "contexulizes" start velocity and potentials.
    // This can be done in the command struct.
    pub fn contextualize(&mut self) -> &mut Self {
        if let Some(atoms) = &mut self.ctx.atoms {
            if atoms.n_atoms == 0 {
                println!("Atoms has not been filled by the input file");
                std::process::exit(1);
            }
            if let Some(starting_velocity) = &self.ctx.starting_velocity {
                if starting_velocity.start_velocity {
                    let start_temperature = match starting_velocity.start_temperature {
                        Some(temperature) => temperature,
                        // If there is no starting temperature given in the input file, a temperature of 300 will be assumed
                        None => 300.0,
                    };
                    let seed = match starting_velocity.seed {
                        Some(seed) => seed,
                        None => 0,
                    };
                    atoms.start_velocities(start_temperature, seed);
                }
            }
        }

        if let Some(potential_args) = &self.ctx.potential_args {
            let mut read_args_style = 0;
            let style = &potential_args.pair_style_args[read_args_style];
            read_args_style += 1;
            match style.as_str() {
                "lj/cut" => {
                    let global_cutoff: f64 = potential_args.pair_style_args[read_args_style]
                        .parse()
                        .unwrap();
                    let mut mgr = LJVOffsetManager::new();

                    for pair_coeff in &potential_args.pair_coeff_args {
                        let mut read_args_coeff = 0;
                        let i: usize = pair_coeff[read_args_coeff].parse().unwrap();
                        read_args_coeff += 1;
                        let j: usize = pair_coeff[read_args_coeff].parse().unwrap();
                        read_args_coeff += 1;
                        let epsilon: f64 = pair_coeff[read_args_coeff].parse().unwrap();
                        read_args_coeff += 1;
                        let sigma: f64 = pair_coeff[read_args_coeff].parse().unwrap();
                        read_args_coeff += 1;
                        let local_rcut: f64 = match pair_coeff.get(read_args_coeff) {
                            Some(rcut) => rcut.parse().unwrap(),
                            None => global_cutoff,
                        };
                        println!("{}, {}, {}", epsilon, sigma, local_rcut);
                        let lj_ij = LennardJones::new(epsilon, sigma, local_rcut, true);
                        mgr.insert((i, j), lj_ij);
                    }
                    self.ctx.mgr = Some(Box::new(mgr));
                }
                _ => print!("Pair Style Unknown"),
            }
        }

        self
    }

    pub fn run(&mut self) {
        self.ctx.run();
    }
}
