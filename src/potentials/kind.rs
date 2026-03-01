use crate::{potentials::{
    lennard_jones::{
        LJManager, LJVOffsetManager, LJVPBuildListManager, LJVParallelManager, LJVerletManager, LennardJones
    },
    potential::{PairPotential, PotentialManager},
}, readers::simulation_context::SimulationContext};

#[allow(dead_code)]
pub enum PotentialManagerKind {
    LJManager(LJManager),
    LJVerletManager(LJVerletManager),
    LJVOffsetManager(LJVOffsetManager),
    LJVParallelManager(LJVParallelManager),
    LJVPBuildListManager(LJVPBuildListManager),
}

impl PotentialManagerKind {
    pub fn run(&self, ctx: &mut SimulationContext) {
        match self {
            Self::LJManager(ljm) => ljm.run(ctx),
            Self::LJVerletManager(ljvm) => ljvm.run(ctx),
            Self::LJVOffsetManager(ljvom) => ljvom.run(ctx),
            Self::LJVParallelManager(ljvpm) => ljvpm.run(ctx),
            Self::LJVPBuildListManager(ljvpbm) => ljvpbm.run(ctx),
        }
    }
}

pub enum PairPotentialKind {
    LennardJones(LennardJones),
}

impl PairPotential for PairPotentialKind {
    fn compute_potential(&self, rij: &na::Vector3<f64>) -> (f64, na::Vector3<f64>) {
        match self {
            PairPotentialKind::LennardJones(lj) => lj.compute_potential(rij),
        }
    }

    fn get_rcut(&self) -> f64 {
        match self {
            PairPotentialKind::LennardJones(lj) => lj.get_rcut(),
        }
    }
}
