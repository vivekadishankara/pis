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
use crate::readers::data_reader::DataReader;

fn main() {
    let args = Args::parse();
    let data_reader = DataReader::new(args.infile);
    let try_compute = false;
    match data_reader.read(args.temperature){
        Ok(mut atoms) => {
            if !try_compute {
                atoms.run(args.timestep, args.steps, &args.dump_path);
            } else {
                let pot_parallel = atoms.compute_potential_verlet_list_parallel();
                println!("{}", pot_parallel);
                let pot_neighbour_list = atoms.compute_potential_neighbour_list_parallel();
                println!("{}", pot_neighbour_list);
            }
        }
        Err(e) => panic!("Could not run the simulation because of {}", e),
    }
}
