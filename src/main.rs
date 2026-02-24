//! Particle Interaction Simulator (PIS) is a rust repository for molecular dynamic(MD)simulation aspiring to be [LAMMPS](<https://www.lammps.org/>) like.
//! 
//! In its current state, PIS implements the pair wise Lennard Jones potential.
//! Because [`nalgebra`] is used to store postions, forces and velocities SIMD acceleration is built into the calculations.
//! 
//! ## Example run
//! To trigger an example calculation run this command
//! ```
//! ./test_command.sh
//! ```
//! An output `dump.lammpstrj` is generated which can be viewed in the OVITO visulation software.
//! 
//! ## Control flow
//! Control in the repository flows through the [`System`] struct. 
//! It is initiated by giving the path to the input file. The system is then read, contextualized and then run to run the simulation.
//! 
//! ```
//! use crate::args_parser::Args;
//! use crate::system::System;
//! 
//! fn main() {
//!     let args = Args::parse();
//!     System::new(args.infile).read().contextualize().run();
//! }
//! ```

extern crate nalgebra as na;

mod args_parser;
mod atoms;
mod constants;
mod ensemble;
mod errors;
mod extensions;
mod math;
mod potentials;
mod readers;
mod simulation_box;
mod system;
mod writers;

use clap::Parser;

use crate::args_parser::Args;
use crate::errors::Result;
use crate::system::System;

fn main() {
    let args = Args::parse();

    if let Err(e) = run(args.infile) {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

fn run(path: String) -> Result<()> {
    System::new(path)
        .read()?
        .contextualize()?
        .run();
    Ok(())
}
