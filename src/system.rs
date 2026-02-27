//! The pivotal struct to initialize and run the system can be found here
use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use crate::{
    errors::{PisError, Result},
    extensions::{ArgsExt, Int32ToUsize},
    potentials::{
        lennard_jones::{LJVOffsetManager, LennardJones},
        potential::PairPotentialManager,
    },
    readers::{input_file::commands::Command, simulation_context::SimulationContext},
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
    ctx: SimulationContext,
}

impl System {
    /// The constructor for the System which takes in the path to the input file
    pub fn new(infile: String) -> Self {
        let ctx = SimulationContext::default();
        Self { infile, ctx }
    }

    /// Reads the input file and collects all the arguments provided by the input file.
    pub fn read(&mut self) -> Result<&mut Self> {
        let file = File::open(&self.infile).map_err(|e| PisError::InputFileError {
            path: self.infile.clone(),
            source: e,
        })?;
        let reader = BufReader::new(file);

        for (line_num, line) in reader.lines().enumerate() {
            let line_num = line_num + 1;
            let line = line.map_err(|e| PisError::DataFileError {
                path: self.infile.clone(),
                line: line_num,
                source: e,
            })?;
            let line = line.trim();

            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            // Considering only the part of the line that is before the commented part.
            let uncommented = line
                .split_once("#")
                .map(|(before, _)| before)
                .unwrap_or(line)
                .trim();

            let line_split: Vec<&str> = uncommented.split_whitespace().collect();

            if line_split.len() < 2 {
                continue;
            }
            let command = line_split[0];
            let args = &line_split[1..];

            match Command::from_str(command) {
                Some(cmd) => cmd.run(args, line_num, &mut self.ctx)?,
                None => {
                    return Err(PisError::UnknownCommand {
                        command: command.to_string(),
                        line: line_num,
                    })
                }
            }
        }
        Ok(self)
    }

    /// Contextualises the system with respect to the arguments read in the [`System::read`].
    /// It actually creates the molecular system to be run.
    // TODO: Maybe the function can be diluted out since it only "contexulizes" start velocity and potentials.
    // This can be done in the command struct.
    pub fn contextualize(&mut self) -> Result<&mut Self> {
        if let Some(atoms) = &mut self.ctx.atoms {
            if atoms.n_atoms == 0 {
                return Err(PisError::NoAtomsDefined);
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
                    atoms.start_velocities(start_temperature, seed)?;
                }
            }
        }

        if let Some(potential_args) = &self.ctx.potential_args {
            let mut read_args_style = 0;
            let style_args: Vec<&str> = potential_args
                .pair_style_args
                .iter()
                .map(|s| s.as_str())
                .collect();
            let style = style_args.get_required(read_args_style, potential_args.pair_style_line)?;
            read_args_style += 1;
            match style {
                "lj/cut" => {
                    let global_cutoff: f64 = style_args
                        .parse_float_at(read_args_style, potential_args.pair_style_line)?;
                    let mut mgr = LJVOffsetManager::new();

                    for (pair_coeff, &coeff_line) in potential_args
                        .pair_coeff_args
                        .iter()
                        .zip(&potential_args.pair_coeff_lines)
                    {
                        let pair_coeff: Vec<&str> = pair_coeff.iter().map(|s| s.as_str()).collect();
                        let mut read_args_coeff = 0;
                        let i: usize = pair_coeff
                            .parse_int_at(read_args_coeff, coeff_line)?
                            .convert_to_usize(coeff_line)?;
                        read_args_coeff += 1;
                        let j: usize = pair_coeff
                            .parse_int_at(read_args_coeff, coeff_line)?
                            .convert_to_usize(coeff_line)?;
                        read_args_coeff += 1;
                        let epsilon: f64 =
                            pair_coeff.parse_float_at(read_args_coeff, coeff_line)?;
                        read_args_coeff += 1;
                        let sigma: f64 = pair_coeff.parse_float_at(read_args_coeff, coeff_line)?;
                        read_args_coeff += 1;
                        let local_rcut: f64 =
                            match pair_coeff.parse_float_at(read_args_coeff, coeff_line) {
                                Ok(rcut) => rcut,
                                Err(_) => global_cutoff,
                            };
                        let lj_ij = LennardJones::new(epsilon, sigma, local_rcut, true);
                        mgr.insert((i, j), lj_ij);
                    }
                    self.ctx.mgr = Some(Box::new(mgr));
                }
                _ => {
                    return Err(PisError::UnknownPairStyle {
                        style: style.to_string(),
                    })
                }
            }
        }

        Ok(self)
    }

    pub fn run(&mut self) {
        self.ctx.run();
    }
}
