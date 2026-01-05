//! Structs to read individual commands in the input file can be found here
use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use na::{DVector, Matrix3, Matrix3xX};

use crate::{
    atoms::new::Atoms,
    potentials::{
        lennard_jones::{LJVOffsetManager, LennardJones},
        potential::PairPotentialManager,
    },
    readers::simulation_context::{
        MTKBarostatArgs, NHThermostatChainArgs, PotentialArgs, SimulationContext, StartVelocity,
        VelocityDistribution,
    },
    simulation_box::SimulationBox,
};

/// The trait which needs to be implemented to write a parser for the arguments of a particular command.
// TODO: It currently returns a anyhow::Result which is not very good at handling errors. This needs to be changed.
// In place, custom errors need to be defined to be returned on encountering different error scenarios
pub trait Command {
    /// The parser function for each command
    fn run(&self, args: &[&str], ctx: &mut SimulationContext) -> anyhow::Result<()>;
}

/// argument parser for the "timestep" command
pub struct TimeStep;

impl Command for TimeStep {
    fn run(&self, args: &[&str], ctx: &mut SimulationContext) -> anyhow::Result<()> {
        ctx.timestep = args[0].parse()?;
        Ok(())
    }
}

/// Argument parser for the "run" command
pub struct RunSteps;

impl Command for RunSteps {
    fn run(&self, args: &[&str], ctx: &mut SimulationContext) -> anyhow::Result<()> {
        ctx.steps = args[0].parse()?;
        Ok(())
    }
}

/// Argument parser for the "velocity" command
pub struct Velocity;

impl Command for Velocity {
    fn run(&self, args: &[&str], ctx: &mut SimulationContext) -> anyhow::Result<()> {
        let mut read_args = 0;
        let mut start_velocity = StartVelocity::default();
        start_velocity.group = String::from(args[read_args]);
        read_args += 1;
        let style = args[read_args];
        read_args += 1;
        match style {
            // Of all the style options in the velocity command only create is available right now
            "create" => {
                start_velocity.start_temperature = Some(args[read_args].parse()?);
                read_args += 1;
                start_velocity.seed = match args.get(read_args) {
                    Some(entry) => {
                        let entry_seed = match entry.parse::<usize>() {
                            Ok(read_seed) => {
                                read_args += 1;
                                read_seed
                            }
                            Err(_) => 0,
                        };
                        Some(entry_seed)
                    }
                    None => Some(0),
                };
            }
            _ => println!("velocity style unknown"),
        }
        loop {
            let keyword = match args.get(read_args) {
                Some(entry) => {
                    read_args += 1;
                    *entry
                }
                None => {
                    break;
                }
            };

            match keyword {
                "dist" => match args[read_args] {
                    "uniform" => start_velocity.dist = Some(VelocityDistribution::Uniform),
                    "gaussian" => start_velocity.dist = Some(VelocityDistribution::Gaussian),
                    _ => println!("velocity distribution unknown"),
                },
                _ => {
                    println!("velocity keyword unknown");
                    break;
                }
            };
        }
        // Currently one one velocity command is allowed in the intput file.
        ctx.starting_velocity = Some(start_velocity);

        Ok(())
    }
}

/// Argument parser for the "read_data" command
pub struct ReadData;

impl Command for ReadData {
    fn run(&self, args: &[&str], ctx: &mut SimulationContext) -> anyhow::Result<()> {
        // the first argument is the path to the data file
        let file = File::open(args[0])?;
        let reader = BufReader::new(file);

        let mut section = String::new();
        let mut n_atoms: usize = 0;

        let (mut xlo, mut xhi): (f64, f64) = (0.0, 1.0);
        let (mut ylo, mut yhi): (f64, f64) = (0.0, 1.0);
        let (mut zlo, mut zhi): (f64, f64) = (0.0, 1.0);

        let mut masses: Vec<f64> = Vec::new();

        let mut type_ids: DVector<usize> = DVector::zeros(0);
        let mut positions: Matrix3xX<f64> = Matrix3xX::zeros(0);
        let mut velocities: Matrix3xX<f64> = Matrix3xX::zeros(0);

        let mut start_velocities = true;

        let mut mgr = LJVOffsetManager::new();

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();

            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let line_split: Vec<&str> = line.split_whitespace().collect();
            // Taking care of the possibility that the parser comes accross the headers in the data file
            match line_split[0] {
                "Masses" | "Atoms" | "PairCoeffs" => {
                    section = line_split[0].to_string();
                    continue;
                }
                "Velocities" => {
                    start_velocities = false;
                    section = line_split[0].to_string();
                    continue;
                }
                _ => {}
            }

            if line_split.len() > 1 {
                match line_split[1] {
                    // reading the number of atoms in the system by reading the number in the line:
                    // "3 atoms"
                    "atoms" => {
                        n_atoms = line_split[0].parse()?;
                        type_ids = DVector::zeros(n_atoms);
                        positions = Matrix3xX::zeros(n_atoms);
                        velocities = Matrix3xX::zeros(n_atoms);
                        continue;
                    }
                    // reading the number of atom types in the system by reading the number in the line:
                    // "2 atom types"
                    "atom" => {
                        let n_types = line_split[0].parse()?;
                        masses.resize(n_types, 0.0);
                        continue;
                    }
                    _ => {}
                }
            }

            if line_split.iter().len() > 2 {
                match line_split[2] {
                    // reading the mlo and mhi numbers (where m can be x, y and z) in the line:
                    // "0.0 5.0 xlo xhi"
                    "xlo" => {
                        xlo = line_split[0].parse()?;
                        xhi = line_split[1].parse()?;
                        continue;
                    }
                    "ylo" => {
                        ylo = line_split[0].parse()?;
                        yhi = line_split[1].parse()?;
                        continue;
                    }
                    "zlo" => {
                        zlo = line_split[0].parse()?;
                        zhi = line_split[1].parse()?;
                        continue;
                    }
                    _ => {}
                };
            }

            match section.as_str() {
                "Masses" => {
                    // reading the masses along with the atoms types
                    let type_id: usize = line_split[0].parse()?;
                    let mass: f64 = line_split[1].parse()?;
                    masses[type_id - 1] = mass;
                }
                "PairCoeffs" => {
                    // reading the pair coefficients for the lennard jones potentials. If the lines looks like this:
                    // "1 0.238 3.405 8.5", this means that the coefficients are for the atoms of type 1 and 1.
                    // if the lines looks like this:
                    // "1 2 0.238 3.405 8.5", then the coefficients are for interations between atom types 1 and 2
                    let i: usize = line_split[0].parse()?;
                    if let Ok(epsilon) = line_split[1].parse::<f64>() {
                        let sigma: f64 = line_split[2].parse()?;
                        let rcut: f64 = match line_split.get(3) {
                            Some(rcut_str) => rcut_str.parse::<f64>()?,
                            None => 2.5 * sigma,
                        };
                        let lj_ii = LennardJones::new(epsilon, sigma, rcut, true);
                        mgr.insert((i, i), lj_ii);
                    } else {
                        let j: usize = line_split[1].parse()?;
                        let epsilon: f64 = line_split[2].parse()?;
                        let sigma: f64 = line_split[3].parse()?;
                        let rcut: f64 = match line_split.get(4) {
                            Some(rcut_str) => rcut_str.parse::<f64>()?,
                            None => 2.5 * sigma,
                        };
                        let lj_ij = LennardJones::new(epsilon, sigma, rcut, true);
                        mgr.insert((i, j), lj_ij);
                    }
                    continue;
                }
                "Atoms" => {
                    // reading atoms from the data file. The line looks like this:
                    // "31 1 0.0 0.0 0.0"
                    // atom_number, atom_type, x_coord, y_coord, z_coord
                    let mut id: usize = line_split[0].parse()?;
                    id -= 1;
                    let type_id: usize = line_split[1].parse()?;
                    type_ids[id] = type_id;
                    let x: f64 = line_split[2].parse()?;
                    let y: f64 = line_split[3].parse()?;
                    let z: f64 = line_split[4].parse()?;
                    positions[(0, id)] = x;
                    positions[(1, id)] = y;
                    positions[(2, id)] = z;
                }
                "Velocities" => {
                    // reading atom velocities from the line that looks like this:
                    // "31 0.0 0.0 0.0"
                    // atom_number, x_coord, y_coord, z_coord
                    let mut id: usize = line_split[0].parse()?;
                    id -= 1;
                    let x: f64 = line_split[1].parse()?;
                    let y: f64 = line_split[2].parse()?;
                    let z: f64 = line_split[3].parse()?;
                    velocities[(0, id)] = x;
                    velocities[(1, id)] = y;
                    velocities[(2, id)] = z;
                }
                _ => {}
            }
        }
        // if atoms already exists in the simulation context
        // TODO: the atoms needs to be added to the original atoms rather tahn replacing them
        if let Some(ctx_atoms) = &mut ctx.atoms {
            ctx_atoms.n_atoms = n_atoms;
            ctx_atoms.type_ids = type_ids;
            ctx_atoms.masses = masses;
            ctx_atoms.positions = positions;
            ctx_atoms.velocities = velocities;
            ctx_atoms.sim_box =
                SimulationBox::from_lammps_data(xlo, xhi, ylo, yhi, zlo, zhi, 0.0, 0.0, 0.0);
        } else {
            // This is the actual place where the atoms read above are assigned to the simulation context
            ctx.atoms = Some(Atoms {
                n_atoms,
                type_ids,
                masses,
                positions,
                velocities,
                forces: Matrix3xX::zeros(n_atoms),
                sim_box: SimulationBox::from_lammps_data(
                    xlo, xhi, ylo, yhi, zlo, zhi, 0.0, 0.0, 0.0,
                ),
            });
        }

        // the potential manager is being assigned here
        if !mgr.is_empty() {
            ctx.mgr = Some(Box::new(mgr));
        }

        // the starting velocites are being assigned here so that they can be initialized in the contextualize part
        if let Some(velocity) = &mut ctx.starting_velocity {
            velocity.start_velocity = start_velocities;
        } else {
            ctx.starting_velocity = Some(StartVelocity {
                group: String::from("all"),
                start_velocity: start_velocities,
                start_temperature: None,
                seed: None,
                dist: None,
            })
        }
        Ok(())
    }
}

/// Argument parser for the "fix" command
pub struct Fix;

impl Command for Fix {
    fn run(&self, args: &[&str], ctx: &mut SimulationContext) -> anyhow::Result<()> {
        let mut read_args: usize = 0;

        let name = String::from(args[read_args]);
        read_args += 1;
        let group = String::from(args[read_args]);
        read_args += 1;
        // Style can either be "npt" or "nvt"
        let style = args[read_args];
        read_args += 1;

        loop {
            let keyword = match args.get(read_args) {
                Some(entry) => {
                    read_args += 1;
                    *entry
                }
                None => {
                    break;
                }
            };
            match keyword {
                "temp" => {
                    let start_temperature: f64 = args[read_args].parse()?;
                    read_args += 1;
                    let end_temperature: f64 = args[read_args].parse()?;
                    read_args += 1;
                    let tau: f64 = args[read_args].parse()?;
                    read_args += 1;
                    let name = name.clone();
                    let group = group.clone();
                    if style == "npt" || style == "nvt" {
                        let nh_chain_args = NHThermostatChainArgs {
                            name,
                            group,
                            start_temperature,
                            end_temperature,
                            tau,
                        };
                        ctx.nh_chain_args = Some(nh_chain_args);
                    }
                }
                "iso" => {
                    let start_pressure: f64 = args[read_args].parse()?;
                    read_args += 1;
                    let _end_pressure: f64 = args[read_args].parse()?;
                    read_args += 1;
                    let tau: f64 = args[read_args].parse()?;
                    read_args += 1;
                    let name = name.clone();
                    let group = group.clone();
                    if style == "npt" {
                        let target_pressure = Matrix3::identity() * start_pressure;
                        let mtk_barostat_args = MTKBarostatArgs {
                            name,
                            group,
                            start_pressure: target_pressure,
                            end_pressure: target_pressure,
                            tau,
                        };
                        ctx.mtk_barostat_args = Some(mtk_barostat_args);
                    }
                }
                _ => println!("Unknow keyword for fix command {}", keyword),
            }
        }
        Ok(())
    }
}


/// Argument parser for the "pair_style" command
pub struct PairStyle;

impl Command for PairStyle {
    fn run(&self, args: &[&str], ctx: &mut SimulationContext) -> anyhow::Result<()> {
        let args: Vec<String> = args.iter().map(|entry| entry.to_string()).collect();
        let mut potential_args = PotentialArgs::default();
        potential_args.pair_style_args = args;
        ctx.potential_args = Some(potential_args);
        Ok(())
    }
}

/// Argument parser for the "pair_coeff" command
pub struct PairCoeff;

impl Command for PairCoeff {
    fn run(&self, args: &[&str], ctx: &mut SimulationContext) -> anyhow::Result<()> {
        let args: Vec<String> = args.iter().map(|entry| entry.to_string()).collect();
        if let Some(potential_args) = &mut ctx.potential_args {
            potential_args.pair_coeff_args.push(args);
        }
        Ok(())
    }
}

/// Argument parser for the "dump" command
pub struct Dump;

impl Command for Dump {
    fn run(&self, args: &[&str], ctx: &mut SimulationContext) -> anyhow::Result<()> {
        let mut read_args = 0;
        ctx.dump_args.name = args[read_args].to_string();
        read_args += 1;
        ctx.dump_args.group = args[read_args].to_string();
        read_args += 1;
        ctx.dump_args.style = args[read_args].to_string();
        read_args += 1;
        ctx.dump_args.dump_step = args[read_args].parse()?;
        read_args += 1;
        ctx.dump_args.file_name = args[read_args].to_string();
        Ok(())
    }
}
