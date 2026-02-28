use crate::potentials::{
    lennard_jones::LennardJones,
    potential::PairPotential,
};

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
