use clap::Parser;

#[derive(Parser)]
#[command(author, version, about)]
pub struct Args {
    #[arg(short, long, default_value_t = String::from("data.txt"))]
    pub infile: String,

    #[arg(short, long, default_value_t = 50)]
    pub steps: usize,

    #[arg(short, long, default_value_t = 0.001)]
    pub timestep: f64,

    #[arg(short = 'T', long, default_value_t = 300.0)]
    pub temperature: f64,

    #[arg(short, long, default_value_t = String::from("dump.lammpstrj"))]
    pub dump_path: String,
}
