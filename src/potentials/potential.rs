use std::collections::HashMap;
use na::Vector3;

pub trait PairPotential: Send + Sync {
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

    pub fn max_rcut(&self) -> f64 {
        let mut max_rcut = 0.0;
        for potential in self.table.values() {
            let rcut = (*potential).get_rcut();
            if max_rcut < rcut { max_rcut = rcut; }
        }
        max_rcut
    }
}
