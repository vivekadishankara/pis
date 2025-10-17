extern crate nalgebra as na;

mod args_parser;
mod readers;
mod atoms;
mod constants;
mod math;
mod simulation_box;
mod potentials;
mod writers;

use clap::Parser;

use crate::args_parser::Args;
use crate::potentials::lennard_jones::{LennardJonesManager, LennardJonesVerletManager, LennardJonesVerletOffsetManager,
    LennardJonesVerletParallelManager};
use crate::potentials::potential::PotentialManager;
use crate::readers::data_reader::DataReader;

fn main() {
    let args = Args::parse();
    let data_reader = DataReader::new(args.infile);
    let potential = "lj_verlet_parallel";
    
    match potential {
        "lj" => {
            match data_reader.read::<LennardJonesManager>(args.temperature){
                Ok((mut atoms, potential_manager)) => {
                    potential_manager.run_nve(&mut atoms, args.timestep, args.steps, &args.dump_path);
                }
                Err(e) => panic!("Could not run the simulation because of {}", e),
            }
        },
        "lj_verlet" => {
            match data_reader.read::<LennardJonesVerletManager>(args.temperature){
                Ok((mut atoms, potential_manager)) => {
                    potential_manager.run_nve(&mut atoms, args.timestep, args.steps, &args.dump_path);
                }
                Err(e) => panic!("Could not run the simulation because of {}", e),
            }
        }
        "lj_verlet_offset" => {
            match data_reader.read::<LennardJonesVerletOffsetManager>(args.temperature){
                Ok((mut atoms, potential_manager)) => {
                    potential_manager.run_nve(&mut atoms, args.timestep, args.steps, &args.dump_path);
                }
                Err(e) => panic!("Could not run the simulation because of {}", e),
            }
        },
        "lj_verlet_parallel" => {
            match data_reader.read::<LennardJonesVerletParallelManager>(args.temperature){
                Ok((mut atoms, potential_manager)) => {
                    potential_manager.run_nve(&mut atoms, args.timestep, args.steps, &args.dump_path);
                }
                Err(e) => panic!("Could not run the simulation because of {}", e),
            }
        },
        _ => panic!("Potential type unknown")
    }
}
