use std::{
    fs::File,
    io::{BufWriter, Result, Write},
};

use crate::{atoms::new::Atoms, simulation_box::SimulationBox};

pub struct DumpTraj {
    out: BufWriter<File>,
}

impl DumpTraj {
    pub fn new(path: &str) -> Result<Self> {
        let file = File::create(path)?;
        Ok(DumpTraj {
            out: BufWriter::new(file),
        })
    }

    pub fn write_timestep(&mut self, step: usize) -> Result<()> {
        writeln!(self.out, "ITEM: TIMESTEP")?;
        writeln!(self.out, "{}", step)?;
        Ok(())
    }

    pub fn write_natoms(&mut self, n_atoms: usize) -> Result<()> {
        writeln!(self.out, "ITEM: NUMBER OF ATOMS")?;
        writeln!(self.out, "{}", n_atoms)?;
        Ok(())
    }

    pub fn write_bounds(&mut self, sim_box: &SimulationBox) -> Result<()> {
        writeln!(self.out, "ITEM: BOX BOUNDS pp pp pp")?;
        writeln!(self.out, "{} {}", 0.0, sim_box.h[(0, 0)])?;
        writeln!(self.out, "{} {}", 0.0, sim_box.h[(1, 1)])?;
        writeln!(self.out, "{} {}", 0.0, sim_box.h[(2, 2)])?;
        Ok(())
    }

    pub fn write_atoms_info(&mut self, atoms: &Atoms) -> Result<()> {
        writeln!(self.out, "ITEM: ATOMS id type x y z")?;
        for i in 0..atoms.n_atoms {
            let position = atoms.positions.column(i);
            writeln!(
                self.out,
                "{} {} {} {} {}",
                i + 1,
                atoms.type_ids[i],
                position[0],
                position[1],
                position[2]
            )?;
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
