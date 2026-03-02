//! Structs for dumping of trajectory output to be found here
use std::{
    fs::File,
    io::{BufWriter, Write},
};

use crate::{
    atoms::new::Atoms, errors::{PisError, Result}, readers::simulation_context::DumpArgs, simulation_box::SimulationBox,
};

/// This struct helps write the trajectory of the simulation as a .traj text file
pub struct DumpTraj {
    out: BufWriter<File>,
    path: String,
}

impl DumpTraj {
    pub fn new(dump_args: &DumpArgs) -> Result<Self> {
        let file = File::create(dump_args.file_name.as_str())
            .map_err(|e| PisError::DumpCreateError { path: dump_args.file_name.clone(), source: e })?;
        Ok(DumpTraj {
            out: BufWriter::new(file),
            path: dump_args.file_name.clone(),
        })
    }

    fn writeln(&mut self, args: std::fmt::Arguments) -> Result<()> {
        writeln!(self.out, "{}", args)
            .map_err(|e| PisError::DumpWriteError { path: self.path.clone(), source: e })
    }

    pub fn write_timestep(&mut self, step: usize) -> Result<()> {
        self.writeln(format_args!("ITEM: TIMESTEP"))?;
        self.writeln(format_args!("{}", step))?;
        Ok(())
    }

    pub fn write_natoms(&mut self, n_atoms: usize) -> Result<()> {
        self.writeln(format_args!("ITEM: NUMBER OF ATOMS"))?;
        self.writeln(format_args!("{}", n_atoms))?;
        Ok(())
    }

    pub fn write_bounds(&mut self, sim_box: &SimulationBox) -> Result<()> {
        self.writeln(format_args!("ITEM: BOX BOUNDS pp pp pp"))?;
        self.writeln(format_args!("{} {}", 0.0, sim_box.h[(0, 0)]))?;
        self.writeln(format_args!("{} {}", 0.0, sim_box.h[(1, 1)]))?;
        self.writeln(format_args!("{} {}", 0.0, sim_box.h[(2, 2)]))?;
        Ok(())
    }

    pub fn write_atoms_info(&mut self, atoms: &Atoms) -> Result<()> {
        self.writeln(format_args!("ITEM: ATOMS id type x y z"))?;
        for i in 0..atoms.n_atoms {
            let position = atoms.positions.column(i);
            self.writeln(format_args!(
                "{} {} {} {} {}",
                i + 1,
                atoms.type_ids[i],
                position[0],
                position[1],
                position[2]
            ))?;
        }
        Ok(())
    }

    pub fn write_step(&mut self, atoms: &Atoms, step: usize) -> Result<()> {
        self.write_timestep(step)?;
        self.write_natoms(atoms.n_atoms)?;
        self.write_bounds(&atoms.sim_box)?;
        self.write_atoms_info(atoms)?;
        Ok(())
    }
}
