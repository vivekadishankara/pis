extern crate nalgebra as na;

mod args_parser;
mod atoms;
mod constants;
mod math;
mod simulation_box;
mod potentials;

use clap::Parser;

use crate::args_parser::Args;

fn main() {
    let args = Args::parse();

    println!("Input file: {}", args.infile);
    println!("Steps: {}", args.steps);
    println!("Timestep: {}", args.timestep);
}
