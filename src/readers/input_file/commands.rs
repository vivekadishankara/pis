use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use na::{DVector, Matrix3, Matrix3xX};

use crate::{
    atoms::new::Atoms, errors::{PisError, Result}, potentials::{
        lennard_jones::{LJVOffsetManager, LennardJones},
        potential::PairPotentialManager,
    }, readers::simulation_context::{
        MTKBarostatArgs, NHThermostatChainArgs, PotentialArgs, SimulationContext, StartVelocity,
        VelocityDistribution,
    }, simulation_box::SimulationBox
};

trait ArgsExt {
    fn get_required(&self, index: usize, line: usize) -> Result<&str>;
    fn parse_int_at(&self, index: usize, line: usize) -> Result<i32>;
    fn parse_float_at(&self, index: usize, line: usize) -> Result<f64>;
}

impl ArgsExt for [&str] {
    fn get_required(&self, index: usize, line: usize) -> Result<&str> {
        self.get(index)
            .copied()
            .ok_or(PisError::MissingArgument { line })
    }
    
    fn parse_int_at(&self, index: usize, line: usize) -> Result<i32> {
        let arg = self.get_required(index, line)?;
        arg.parse()
            .map_err(|e| PisError::IntParseError {
                string: arg.to_string(),
                source: e,
            })
    }

    fn parse_float_at(&self, index: usize, line: usize) -> Result<f64> {
        let arg = self.get_required(index, line)?;
        arg.parse()
            .map_err(|e| PisError::FloatParseError { 
                string: arg.to_string(),
                source: e 
            })
    }
}

trait Int32ToUsize {
    fn convert_to_usize(&self, line: usize) -> Result<usize>;
}

impl Int32ToUsize for i32 {
    fn convert_to_usize(&self, line: usize) -> Result<usize> {
        (*self).try_into()
            .map_err(|_| PisError::NegativeValue{
                value: (*self),
                line: line
            })
    }
}

pub trait Command {
    fn run(&self, args: &[&str], line:usize, ctx: &mut SimulationContext) -> Result<()>;
}

pub struct TimeStep;

impl Command for TimeStep {
    fn run(&self, args: &[&str], line:usize, ctx: &mut SimulationContext) -> Result<()> {
        ctx.timestep = args.parse_float_at(0, line)?;
        Ok(())
    }
}

pub struct RunSteps;

impl Command for RunSteps {
    fn run(&self, args: &[&str], line:usize, ctx: &mut SimulationContext) -> Result<()> {
        ctx.steps = args.parse_int_at(0, line)?.convert_to_usize(line)?;
        Ok(())
    }
}

pub struct Velocity;

impl Command for Velocity {
    fn run(&self, args: &[&str], line:usize, ctx: &mut SimulationContext) -> Result<()> {
        let mut read_args = 0;
        let mut start_velocity = StartVelocity::default();
        start_velocity.group = String::from(args.get_required(read_args, line)?);
        read_args += 1;
        let style = args.get_required(read_args, line)?;
        read_args += 1;
        match style {
            "create" => {
                start_velocity.start_temperature = Some(args.parse_float_at(read_args, line)?);
                read_args += 1;
                start_velocity.seed = match args.parse_int_at(read_args, line) {
                    Ok(seed) => {
                        read_args += 1;
                        Some(seed.convert_to_usize(line)?)
                    },
                    Err(_) => Some(0),
                };
            }
            _ => return Err(PisError::InvalidArgument { string: style.to_string(), line }),
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
                "dist" => {
                    let keyword_arg = args.get_required(read_args, line)?;
                    match keyword_arg {
                    "uniform" => start_velocity.dist = Some(VelocityDistribution::Uniform),
                    "gaussian" => start_velocity.dist = Some(VelocityDistribution::Gaussian),
                    _ => return Err(PisError::InvalidArgument { string: keyword_arg.to_string(), line }),
                    }
                },
                
                _ => return Err(PisError::InvalidArgument { string: keyword.to_string(), line })
            };
        }
        ctx.starting_velocity = Some(start_velocity);

        Ok(())
    }
}

pub struct ReadData;

impl Command for ReadData {
    fn run(&self, args: &[&str], line:usize, ctx: &mut SimulationContext) -> Result<()> {
        let path = args.get_required(0, line)?;
        let file = File::open(path)
            .map_err(|e| PisError::InputFileError { path: path.to_string(), source: e })?;
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

        for (line_num, line) in reader.lines().enumerate() {
            let line = line
                .map_err(|e| PisError::DataFileError { 
                    path: path.to_string(), 
                    line: line_num, 
                    source: e 
                })?;
            let line = line.trim();

            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let line_split: Vec<&str> = line.split_whitespace().collect();
            
            let first_word = line_split.get_required(0, line_num)?;
            match first_word {
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
                    "atoms" => {
                        n_atoms = line_split.parse_int_at(0, line_num)?.convert_to_usize(line_num)?;
                        type_ids = DVector::zeros(n_atoms);
                        positions = Matrix3xX::zeros(n_atoms);
                        velocities = Matrix3xX::zeros(n_atoms);
                        continue;
                    }
                    "atom" => {
                        let n_types = line_split.parse_int_at(0, line_num)?.convert_to_usize(line_num)?;
                        masses.resize(n_types, 0.0);
                        continue;
                    }
                    _ => {}
                }
            }

            if line_split.iter().len() > 2 {
                match line_split[2] {
                    "xlo" => {
                        xlo = line_split.parse_float_at(0, line_num)?;
                        xhi = line_split.parse_float_at(1, line_num)?;
                        continue;
                    }
                    "ylo" => {
                        ylo = line_split.parse_float_at(0, line_num)?;
                        yhi = line_split.parse_float_at(1, line_num)?;
                        continue;
                    }
                    "zlo" => {
                        zlo = line_split.parse_float_at(0, line_num)?;
                        zhi = line_split.parse_float_at(1, line_num)?;
                        continue;
                    }
                    _ => {}
                };
            }

            match section.as_str() {
                "Masses" => {
                    let type_id: usize = line_split.parse_int_at(0, line_num)?.convert_to_usize(line_num)?;
                    let mass: f64 = line_split.parse_float_at(1, line_num)?;
                    if type_id < 1 {
                        return Err(PisError::InvalidAtomType { type_id })
                    }
                    masses[type_id - 1] = mass;
                }
                "PairCoeffs" => {
                    let i: usize = line_split.parse_int_at(0, line_num)?.convert_to_usize(line_num)?;
                    if let Ok(epsilon) = line_split.parse_float_at(1, line_num) {
                        let sigma: f64 = line_split.parse_float_at(2, line_num)?;
                        let rcut: f64 = match line_split.parse_float_at(3, line_num) {
                            Ok(rcut) => rcut,
                            Err(_) => 2.5 * sigma,
                        };
                        let lj_ii = LennardJones::new(epsilon, sigma, rcut, true);
                        mgr.insert((i, i), lj_ii);
                    } else {
                        let j: usize = line_split.parse_int_at(1, line_num)?.convert_to_usize(line_num)?;
                        let epsilon: f64 = line_split.parse_float_at(2, line_num)?;
                        let sigma: f64 = line_split.parse_float_at(3, line_num)?;
                        let rcut: f64 = match line_split.parse_float_at(4, line_num) {
                            Ok(rcut) => rcut,
                            Err(_) => 2.5 * sigma,
                        };
                        let lj_ij = LennardJones::new(epsilon, sigma, rcut, true);
                        mgr.insert((i, j), lj_ij);
                    }
                    // let mut j: usize;
                    continue;
                }
                "Atoms" => {
                    let mut id: usize = line_split.parse_int_at(0, line_num)?.convert_to_usize(line_num)?;
                    if id == 0 || id > n_atoms {
                        return Err(PisError::AtomCountMismatch { expected: n_atoms, found: id })
                    }
                    id -= 1;
                    let type_id: usize = line_split.parse_int_at(1, line_num)?.convert_to_usize(line_num)?;
                    type_ids[id] = type_id;
                    let x: f64 = line_split.parse_float_at(2, line_num)?;
                    let y: f64 = line_split.parse_float_at(3, line_num)?;
                    let z: f64 = line_split.parse_float_at(4, line_num)?;
                    positions[(0, id)] = x;
                    positions[(1, id)] = y;
                    positions[(2, id)] = z;
                }
                "Velocities" => {
                    let mut id: usize = line_split.parse_int_at(0, line_num)?.convert_to_usize(line_num)?;
                    if id == 0 || id > n_atoms {
                        return Err(PisError::AtomCountMismatch { expected: n_atoms, found: id })
                    }
                    id -= 1;
                    let x: f64 = line_split.parse_float_at(1, line_num)?;
                    let y: f64 = line_split.parse_float_at(2, line_num)?;
                    let z: f64 = line_split.parse_float_at(3, line_num)?;
                    velocities[(0, id)] = x;
                    velocities[(1, id)] = y;
                    velocities[(2, id)] = z;
                }
                _ => {}
            }
        }

        if let Some(ctx_atoms) = &mut ctx.atoms {
            ctx_atoms.n_atoms = n_atoms;
            ctx_atoms.type_ids = type_ids;
            ctx_atoms.masses = masses;
            ctx_atoms.positions = positions;
            ctx_atoms.velocities = velocities;
            ctx_atoms.sim_box =
                SimulationBox::from_lammps_data(xlo, xhi, ylo, yhi, zlo, zhi, 0.0, 0.0, 0.0);
        } else {
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

        if !mgr.is_empty() {
            ctx.mgr = Some(Box::new(mgr));
        }

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

pub struct Fix;

impl Command for Fix {
    fn run(&self, args: &[&str], line:usize, ctx: &mut SimulationContext) -> Result<()> {
        let mut read_args: usize = 0;

        let name = String::from(args.get_required(read_args, line)?);
        read_args += 1;
        let group = String::from(args.get_required(read_args, line)?);
        read_args += 1;
        let style = args.get_required(read_args, line)?;
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
                    let start_temperature: f64 = args.parse_float_at(read_args, line)?;
                    read_args += 1;
                    let end_temperature: f64 = args.parse_float_at(read_args, line)?;
                    read_args += 1;
                    let tau: f64 = args.parse_float_at(read_args, line)?;
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
                    let start_pressure: f64 = args.parse_float_at(read_args, line)?;
                    read_args += 1;
                    let _end_pressure: f64 = args.parse_float_at(read_args, line)?;
                    read_args += 1;
                    let tau: f64 = args.parse_float_at(read_args, line)?;
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

pub struct PairStyle;

impl Command for PairStyle {
    fn run(&self, args: &[&str], _line:usize, ctx: &mut SimulationContext) -> Result<()> {
        let args: Vec<String> = args.iter().map(|entry| entry.to_string()).collect();
        let mut potential_args = PotentialArgs::default();
        potential_args.pair_style_args = args;
        ctx.potential_args = Some(potential_args);
        Ok(())
    }
}

pub struct PairCoeff;

impl Command for PairCoeff {
    fn run(&self, args: &[&str], _line:usize, ctx: &mut SimulationContext) -> Result<()> {
        let args: Vec<String> = args.iter().map(|entry| entry.to_string()).collect();
        if let Some(potential_args) = &mut ctx.potential_args {
            potential_args.pair_coeff_args.push(args);
        }
        Ok(())
    }
}

pub struct Dump;

impl Command for Dump {
    fn run(&self, args: &[&str], line:usize, ctx: &mut SimulationContext) -> Result<()> {
        let mut read_args = 0;
        ctx.dump_args.name = args.get_required(read_args, line)?.to_string();
        read_args += 1;
        ctx.dump_args.group = args.get_required(read_args, line)?.to_string();
        read_args += 1;
        ctx.dump_args.style = args.get_required(read_args, line)?.to_string();
        read_args += 1;
        ctx.dump_args.dump_step = args.parse_int_at(read_args, line)?.convert_to_usize(line)?;
        read_args += 1;
        ctx.dump_args.file_name = args.get_required(read_args, line)?.to_string();
        Ok(())
    }
}
