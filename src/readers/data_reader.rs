use std::{fs::File, io::{BufRead, BufReader}, usize};

use na::{DVector, Matrix3xX};

use crate::{atoms::new::Atoms, potentials::{lennard_jones::{LennardJones}, potential::PairPotentialManager}, simulation_box::SimulationBox};

pub struct DataReader {
    infile: String,
}

impl DataReader {
    pub fn new(infile: String) -> Self {
        Self { infile }
    }

    pub fn read<T: PairPotentialManager>(&self, temperature: f64) -> anyhow::Result<(Atoms, T)> {
        let file = File::open(&self.infile)?;
        let reader = BufReader::new(file);

        let mut section = String::new();
        let mut n_atoms: usize = 0;

        let (mut xlo, mut xhi): (f64, f64) = (0.0, 1.0);
        let (mut ylo, mut yhi): (f64, f64) = (0.0, 1.0);
        let (mut zlo, mut zhi): (f64, f64) = (0.0, 1.0);

        let mut masses: Vec<f64> = Vec::new();

        let mut type_ids: DVector<usize> = DVector::zeros(0);
        let mut positions: Matrix3xX<f64> = Matrix3xX::zeros(0);
        let mut velocities: Matrix3xX<f64> = Matrix3xX::zeros(0);

        let mut start_velocities = true;

        let mut mgr = T::new();

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();

            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let line_split: Vec<&str> = line.split_whitespace().collect();

            match line_split[0] {
                "Masses" | "Atoms" | "PairCoeffs" => {
                    section = line_split[0].to_string();
                    continue;
                },
                "Velocities" => {
                    start_velocities = false;
                    section = line_split[0].to_string();
                    continue;
                }
                _ => {}
            }

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
                        let n_types = line_split[0].parse()?;
                        masses.resize(n_types, 0.0);
                        continue;
                    },
                    _ => {},
                }
            }

            if line_split.iter().len() > 2 {
                match line_split[2] {
                    "xlo" => {
                        xlo = line_split[0].parse()?;
                        xhi = line_split[1].parse()?;
                        continue;
                    },
                    "ylo" => {
                        ylo = line_split[0].parse()?;
                        yhi = line_split[1].parse()?;
                        continue;
                    },
                    "zlo" => {
                        zlo = line_split[0].parse()?;
                        zhi = line_split[1].parse()?;
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
                    if let Ok(epsilon) = line_split[1].parse::<f64>() {
                        let sigma: f64 = line_split[2].parse()?;
                        let rcut: f64 = match line_split.get(3) {
                            Some(rcut_str) => rcut_str.parse::<f64>()?,
                            None => 2.5 * sigma,
                        };
                        let lj_ii = LennardJones::new(epsilon, sigma, rcut, true);
                        mgr.insert((i,i), lj_ii);
                    } else {
                        let j: usize = line_split[1].parse()?;
                        let epsilon: f64 = line_split[2].parse()?;
                        let sigma: f64 = line_split[3].parse()?;
                        let rcut: f64 = match line_split.get(4) {
                            Some(rcut_str) => rcut_str.parse::<f64>()?,
                            None => 2.5 * sigma,
                        };
                        let lj_ij = LennardJones::new(epsilon, sigma, rcut, true);
                        mgr.insert((i,j), lj_ij);
                    }
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
        // println!("{} {} {:?} {:?} {:?} {} {} {}", n_atoms, n_types, x_bounds, y_bounds, z_bounds, type_ids, positions, velocities);
        let mut atoms = Atoms{
            n_atoms,
            type_ids,
            masses,
            positions,
            velocities,
            forces: Matrix3xX::zeros(n_atoms),
            sim_box: SimulationBox::from_lammps_data(xlo, xhi, ylo, yhi, zlo, zhi, 0.0, 0.0, 0.0),
        };

        if start_velocities {
            atoms.start_velocities(temperature);
        }
        Ok((atoms, mgr))
    }
    
}
