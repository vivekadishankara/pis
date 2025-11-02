use clap::Parser;

#[derive(Parser)]
#[command(author, version, about)]
pub struct Args {
    #[arg(short, long, default_value_t = String::from("input.pis"))]
    pub infile: String,
}
