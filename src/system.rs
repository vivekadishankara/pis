use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
};

use crate::readers::{
    input_file::commands::{Command, Fix, ReadData, RunSteps, TimeStep, Velocity},
    simulation_context::SimulationContext,
};

pub struct System {
    infile: String,
    command_hash: HashMap<String, Box<dyn Command>>,
    ctx: SimulationContext,
}

impl System {
    pub fn new(infile: String) -> Self {
        let command_hash = HashMap::new();
        let ctx = SimulationContext::default();
        Self {
            infile,
            command_hash,
            ctx,
        }
    }

    fn build_command_hash(&mut self) {
        self.command_hash
            .insert(String::from("timestep"), Box::new(TimeStep));
        self.command_hash
            .insert(String::from("run"), Box::new(RunSteps));
        self.command_hash
            .insert(String::from("velocity"), Box::new(Velocity));
        self.command_hash
            .insert(String::from("read_data"), Box::new(ReadData));
        self.command_hash.insert(String::from("fix"), Box::new(Fix));
    }

    pub fn read(&mut self) -> &mut Self {
        let file = File::open(&self.infile).expect("Failed to open {infile}");
        let reader = BufReader::new(file);

        self.build_command_hash();

        for line in reader.lines() {
            let line = line.unwrap();
            let line = line.trim();

            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let uncommented = line.splitn(2, '#').next().unwrap().trim();

            let line_split: Vec<&str> = uncommented.split_whitespace().collect();

            let command = line_split[0];
            let args = &line_split[1..];

            if let Some(handler) = self.command_hash.get(command) {
                let _ = handler.run(args, &mut self.ctx);
            } else {
                println!("Command Unkown");
            }
        }
        self
    }

    pub fn contextualize(&mut self) -> &mut Self {
        if let Some(atoms) = &mut self.ctx.atoms {
            if atoms.n_atoms == 0 {
                panic!("Atoms has not been filled by the input file")
            }
            if let Some(starting_velocity) = &self.ctx.starting_velocity {
                if starting_velocity.start_velocity {
                    let start_temperature = match starting_velocity.start_temperature {
                        Some(temperature) => temperature,
                        None => 300.0,
                    };
                    let seed = match starting_velocity.seed {
                        Some(seed) => seed,
                        None => 0,
                    };
                    atoms.start_velocities(start_temperature, seed);
                }
            }
        }

        self
    }

    pub fn run(&mut self) {
        self.ctx.run();
    }
}
