use clap::Parser;

/// The command line argument parser for the molecular dynamics run derived from the clap Parser
#[derive(Parser)]
#[command(author, version, about)]
pub struct Args {
    #[arg(short, long, default_value_t = String::from("input.pis"))]
    pub infile: String,
}
