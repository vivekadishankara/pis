use std::collections::HashMap;
use na::Vector3;

pub trait PairPotential {
    fn compute_potetial(&self, rij: &Vector3<f64>) -> (f64, Vector3<f64>);
    fn get_rcut(&self) -> f64;
}

type AtomPair = (usize, usize);

pub struct PairPotentialManager {
    table: HashMap<AtomPair, Box<dyn PairPotential>>,
}

impl PairPotentialManager {
    pub fn new() -> Self {
        Self { table: HashMap::new() }
    }

    pub fn insert<P>(&mut self, key: AtomPair, potential: P) 
    where 
        P: PairPotential + 'static,
    {
        self.table.insert(key, Box::new(potential));
    }

    pub fn get(&self, key: &AtomPair) -> Option<&dyn PairPotential> {
        self.table.get(key).map(|b| b.as_ref())
    }
}
