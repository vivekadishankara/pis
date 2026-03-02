use crate::{
    atoms::new::Atoms, ensemble::{npt::MTKBarostat, nvt::NHThermostatChain}, errors::{PisError, Result}, potentials::{kind::PotentialManagerKind, potential::PotentialManager}, readers::simulation_context::SimulationContext, writers::dump_traj::DumpTraj
};

pub struct Simulation;

impl Simulation {
    pub fn run(
        mgr: &PotentialManagerKind,
        ctx: &mut SimulationContext,
    ) -> Result<()> {
        let ensemble = Ensemble::from_ctx(ctx);
        let dt = ctx.timestep;
        let steps = ctx.steps;
        let dump_step = ctx.dump_args.dump_step;

        let mut nose_hoover_chain = NHThermostatChain::new_from_args(&ctx.nh_chain_args);
        let mut mtk_barostat = MTKBarostat::new_from_args(
            &ctx.mtk_barostat_args,
            &ctx.nh_chain_args,
            ctx.atoms.as_ref().ok_or(PisError::NoAtomsDefined)?.n_atoms
        );

        let mut dumper = DumpTraj::new(&ctx.dump_args)?;
        let atoms = ctx.atoms.as_mut().ok_or(PisError::NoAtomsDefined)?;

        dumper.write_step(atoms, 0)?;
        let first_potential = mgr.compute_potential(atoms);
        println!("{} {}", 0, first_potential);

        for i in 0..steps {
            let atoms = ctx.atoms.as_mut().ok_or(PisError::NoAtomsDefined)?;
            let step_potential = Self::step(mgr, atoms, ensemble, dt, &mut nose_hoover_chain, &mut mtk_barostat, i, steps);

            let atoms = ctx.atoms.as_mut().ok_or(PisError::NoAtomsDefined)?;
            Self::output(atoms, &mut dumper, ensemble, &nose_hoover_chain, &mtk_barostat, i, step_potential, dump_step)?;
        }
        Ok(())
    }

    fn step(
        mgr: &PotentialManagerKind,
        atoms: &mut Atoms,
        ensemble: Ensemble,
        dt: f64,
        nhc: &mut Option<NHThermostatChain>,
        mtk: &mut Option<MTKBarostat>,
        i: usize,
        steps: usize,
    ) -> f64 {
        match ensemble {
            Ensemble::NVE => mgr.verlet_step_nve(atoms, dt),
            Ensemble::NVT => {
                let potential = mgr.verlet_step_nvt_nhc(atoms, dt, nhc.as_mut().unwrap());
                nhc.as_mut().unwrap().calculate_target_temperature(i, steps);
                potential
            }
            Ensemble::NPT => {
                let potential = mgr.verlet_step_npt_mtk(atoms, dt, mtk.as_mut().unwrap(), nhc.as_mut().unwrap());
                nhc.as_mut().unwrap().calculate_target_temperature(i, steps);
                potential
            }
        }
    }

    fn output(
        atoms: &Atoms,
        dumper: &mut DumpTraj,
        ensemble: Ensemble,
        nhc: &Option<NHThermostatChain>,
        mtk: &Option<MTKBarostat>,
        i: usize,
        step_potential: f64,
        dump_step: usize,
    ) -> Result<()> {
        if ((i + 1) % dump_step) == 0 {
            dumper.write_step(atoms, i + 1)?;
        }
        let kinetic_energy = atoms.kinetic_energy();
        let hamiltonian = Self::compute_hamiltonian(ensemble, step_potential, kinetic_energy, nhc, mtk, atoms);
        let temperature = atoms.temerature(kinetic_energy);
        let pressure = atoms.pressure(kinetic_energy);
        println!(
            "{} {:.3} {:.3} {:.3} {:.3} {:.3}",
            i + 1, step_potential, kinetic_energy, hamiltonian, temperature, pressure
        );
        Ok(())
    }

    fn compute_hamiltonian(
        ensemble: Ensemble,
        potential: f64,
        kinetic: f64,
        nhc: &Option<NHThermostatChain>,
        mtk: &Option<MTKBarostat>,
        atoms: &Atoms,
    ) -> f64 {
        let basic = potential + kinetic;
        match ensemble {
            Ensemble::NVE => basic,
            Ensemble::NVT => {
                let nhc = nhc.as_ref().unwrap();
                basic + nhc.kinetic_energy() + nhc.potential_energy(atoms.n_atoms)
            }
            Ensemble::NPT => {
                let nhc = nhc.as_ref().unwrap();
                let mtk = mtk.as_ref().unwrap();
                basic
                    + nhc.kinetic_energy()
                    + nhc.potential_energy(atoms.n_atoms)
                    + mtk.kinetic_energy()
                    + mtk.potential_energy(&atoms.sim_box.h)
            }
        }
    }
}

#[derive(Clone, Copy)]
pub enum Ensemble {
    NVE,
    NVT,
    NPT,
}

impl Ensemble {
    fn from_ctx(ctx: &SimulationContext) -> Self {
        match (&ctx.nh_chain_args, &ctx.mtk_barostat_args) {
            (Some(_), Some(_)) => Ensemble::NPT,
            (Some(_), None) => Ensemble::NVT,
            _ => Ensemble::NVE,
        }
    }
}
