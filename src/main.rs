extern crate nalgebra as na;

mod args_parser;
mod atoms;
mod constants;
mod ensemble;
mod errors;
mod math;
mod potentials;
mod readers;
mod simulation_box;
mod system;
mod writers;

use clap::Parser;

use crate::args_parser::Args;
use crate::system::System;

fn main() {
    let args = Args::parse();
    System::new(args.infile).read().contextualize().run();
}
