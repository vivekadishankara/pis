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
    match data_reader.read(args.temperature){
        Ok(mut atoms) => {
            atoms.run(args.timestep, args.steps, &args.dump_path);
        }
        Err(e) => panic!("Could not run the simulation because of {}", e),
    }
}
