extern crate nalgebra as na;

mod args_parser;
mod readers;
mod atoms;
mod constants;
mod math;
mod simulation_box;
mod potentials;

use clap::Parser;

use crate::args_parser::Args;
use crate::readers::data_reader::DataReader;

fn main() {
    let args = Args::parse();
    let data_reader = DataReader::new(args.infile);
    let _ = data_reader.read();
}
