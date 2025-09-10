use std::{fs::File, io::{BufRead, BufReader}, usize};

use na::{DVector, Matrix3xX};

pub struct DataReader {
    infile: String,
}

impl DataReader {
    pub fn new(infile: String) -> Self {
        Self { infile }
    }

    pub fn read(&self) -> anyhow::Result<()> {
        let file = File::open(&self.infile)?;
        let reader = BufReader::new(file);

        let mut section = String::new();
        let mut n_atoms: usize = 0;
        let mut n_types: usize = 0;

        let mut x_bounds: Vec<f64> = Vec::with_capacity(2);
        let mut y_bounds: Vec<f64> = Vec::with_capacity(2);
        let mut z_bounds: Vec<f64> = Vec::with_capacity(2);

        let mut masses: Vec<f64> = Vec::new();

        let mut type_ids: DVector<usize> = DVector::zeros(0);
        let mut positions: Matrix3xX<f64> = Matrix3xX::zeros(0);
        let mut velocities: Matrix3xX<f64> = Matrix3xX::zeros(0);

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();

            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let line_split: Vec<&str> = line.split_whitespace().collect();

            match line_split[0] {
                "Masses" | "Atoms" | "Velocities" | "PairCoeffs" => {
                    section = line.to_string();
                    continue;
                }
                _ => {}
            }
            println!("{}", section);

            if line_split.len() > 1 {
                match line_split[1] {
                    "atoms" => {
                        n_atoms = line_split[0].parse()?;
                        type_ids = DVector::zeros(n_atoms);
                        positions= Matrix3xX::zeros(n_atoms);
                        velocities= Matrix3xX::zeros(n_atoms);
                        continue;
                    },
                    "atom" => {
                        n_types = line_split[0].parse()?;
                        masses.resize(n_types, 0.0);
                        continue;
                    },
                    _ => {},
                }
            }

            if line_split.iter().len() > 2 {
                match line_split[2] {
                    "xlo" => {
                        x_bounds = line_split[0..2]
                            .iter()
                            .map(|value| value.parse::<f64>())
                            .collect::<Result<Vec<_>, _>>()?;
                        continue;
                    },
                    "ylo" => {
                        y_bounds = line_split[0..2]
                            .iter()
                            .map(|value| value.parse::<f64>())
                            .collect::<Result<Vec<_>, _>>()?;
                        continue;
                    },
                    "zlo" => {
                        z_bounds = line_split[0..2]
                            .iter()
                            .map(|value| value.parse::<f64>())
                            .collect::<Result<Vec<_>, _>>()?;
                        continue;
                    },
                    _ => {},
                };
            }

            match section.as_str() {
                "Masses" => {
                    let type_id: usize = line_split[0].parse()?;
                    let mass: f64 = line_split[1].parse()?;
                    masses[type_id - 1] = mass;
                },
                "PairCoeffs" => {
                    let i: usize = line_split[0].parse()?;
                    // let mut j: usize;
                    continue;
                },
                "Atoms" => {
                    let mut id: usize = line_split[0].parse()?;
                    id -= 1;
                    let type_id: usize = line_split[1].parse()?;
                    type_ids[id] = type_id;
                    let x: f64 = line_split[2].parse()?;
                    let y: f64 = line_split[3].parse()?;
                    let z: f64 = line_split[4].parse()?;
                    positions[(0, id)] = x;
                    positions[(1, id)] = y;
                    positions[(2, id)] = z;
                },
                "Velocities" => {
                    let mut id: usize = line_split[0].parse()?;
                    id -= 1;
                    let x: f64 = line_split[1].parse()?;
                    let y: f64 = line_split[2].parse()?;
                    let z: f64 = line_split[3].parse()?;
                    velocities[(0, id)] = x;
                    velocities[(1, id)] = y;
                    velocities[(2, id)] = z;
                },
                _ => {},
            }
        }
        println!("{} {} {:?} {:?} {:?} {} {} {}", n_atoms, n_types, x_bounds, y_bounds, z_bounds, type_ids, positions, velocities);

        Ok(())
    }
    
}
