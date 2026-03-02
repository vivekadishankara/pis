use crate::{
    atoms::new::Atoms,
    potentials::{
        lennard_jones::{
            LJManager, LJVOffsetManager, LJVPBuildListManager, LJVParallelManager, LJVerletManager,
            LennardJones,
        },
        potential::{PairPotential, PotentialManager},
    },
};

#[allow(dead_code)]
pub enum PotentialManagerKind {
    LJManager(LJManager),
    LJVerletManager(LJVerletManager),
    LJVOffsetManager(LJVOffsetManager),
    LJVParallelManager(LJVParallelManager),
    LJVPBuildListManager(LJVPBuildListManager),
}

impl PotentialManager for PotentialManagerKind {
    fn compute_potential(&self, atoms: &mut Atoms) -> f64 {
        match self {
            Self::LJManager(ljm) => ljm.compute_potential(atoms),
            Self::LJVerletManager(ljvm) => ljvm.compute_potential(atoms),
            Self::LJVOffsetManager(ljvom) => ljvom.compute_potential(atoms),
            Self::LJVParallelManager(ljvpm) => ljvpm.compute_potential(atoms),
            Self::LJVPBuildListManager(ljvpbm) => ljvpbm.compute_potential(atoms),
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
