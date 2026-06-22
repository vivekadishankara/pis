mod args_parser;

use clap::Parser;
use pis::errors::Result;
use pis::system::System;

use crate::args_parser::Args;

fn main() {
    let args = Args::parse();

    if let Err(e) = run(args.infile) {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

fn run(path: String) -> Result<()> {
    System::new(path).read()?.contextualize()?.run()?;
    Ok(())
}
