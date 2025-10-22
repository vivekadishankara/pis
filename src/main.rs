extern crate nalgebra as na;

mod args_parser;
mod atoms;
mod constants;
mod ensemble;
mod math;
mod potentials;
mod readers;
mod simulation_box;
mod writers;

use clap::Parser;

use crate::args_parser::Args;
use crate::potentials::lennard_jones::{
    LJManager, LJVOffsetManager, LJVPBuildListManager, LJVParallelManager, LJVerletManager,
};
use crate::potentials::potential::PotentialManager;
use crate::readers::data_reader::DataReader;

fn main() {
    let args = Args::parse();
    let data_reader = DataReader::new(args.infile);
    let potential = "lj_verlet_offset";

    match potential {
        "lj" => match data_reader.read::<LJManager>(args.temperature) {
            Ok((mut atoms, potential_manager)) => {
                potential_manager.run(&mut atoms, args.timestep, args.steps, &args.dump_path);
            }
            Err(e) => panic!("Could not run the simulation because of {}", e),
        },
        "lj_verlet" => match data_reader.read::<LJVerletManager>(args.temperature) {
            Ok((mut atoms, potential_manager)) => {
                potential_manager.run(&mut atoms, args.timestep, args.steps, &args.dump_path);
            }
            Err(e) => panic!("Could not run the simulation because of {}", e),
        },
        "lj_verlet_offset" => match data_reader.read::<LJVOffsetManager>(args.temperature) {
            Ok((mut atoms, potential_manager)) => {
                potential_manager.run(&mut atoms, args.timestep, args.steps, &args.dump_path);
            }
            Err(e) => panic!("Could not run the simulation because of {}", e),
        },
        "lj_verlet_parallel" => match data_reader.read::<LJVParallelManager>(args.temperature) {
            Ok((mut atoms, potential_manager)) => {
                potential_manager.run(&mut atoms, args.timestep, args.steps, &args.dump_path);
            }
            Err(e) => panic!("Could not run the simulation because of {}", e),
        },
        "lj_build_list_parallel" => {
            match data_reader.read::<LJVPBuildListManager>(args.temperature) {
                Ok((mut atoms, potential_manager)) => {
                    potential_manager.run(&mut atoms, args.timestep, args.steps, &args.dump_path);
                }
                Err(e) => panic!("Could not run the simulation because of {}", e),
            }
        }
        _ => panic!("Potential type unknown"),
    }
}
