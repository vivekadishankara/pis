//! The member functions for Atoms which relate to neighbour list formation are defined here
#![allow(dead_code)]
use crate::atoms::new::Atoms;

impl Atoms {
    /// This function divides the atoms into cells according to the number of cells in each direction.
    /// The position of atom is first converted into unit box coordinates. Then they are assigned to a cell
    /// according to the cell number computed by using the standard formula
    pub fn rcut_cells(&self, nx: usize, ny: usize, nz: usize) -> Vec<Vec<usize>> {
        let ncell_total = nx * ny * nz;

        let approx_atoms_per_cell = self.n_atoms / ncell_total + 1;
        // The neigbour list is initiated as a vector of vectors. The inner vector is initialized with capacity of atoms per box.
        let mut cells: Vec<Vec<usize>> =
            vec![Vec::with_capacity(approx_atoms_per_cell); ncell_total];

        for (i, r_i) in self.positions.column_iter().enumerate() {
            // The absolute coordinate of the atom is coverted into unit box one
            let s_i = self.sim_box.h_inv * r_i;
            // Each coordinate is then multiplied with the number of cells in that direction
            let mut cx_i = (s_i[0] * nx as f64).floor() as usize;
            let mut cy_i = (s_i[1] * ny as f64).floor() as usize;
            let mut cz_i = (s_i[2] * nz as f64).floor() as usize;

            // Application of periodic boundary conditions as per the simulation box
            if self.sim_box.pbc[0] {
                cx_i = cx_i % nx;
            }
            if self.sim_box.pbc[1] {
                cy_i = cy_i % ny;
            }
            if self.sim_box.pbc[2] {
                cz_i = cz_i % nz;
            }

            // Then the cell numbers cx, cy znd cz are combined into a number according to a standard formula
            let cell_index = Self::cell_index(cx_i, cy_i, cz_i, nx, ny);
            cells[cell_index].push(i);
        }

        cells
    }

    /// This functions takes in an cutoff radius and divides the box into cells according to that radius.
    pub fn divide_into_cells(&self, rcut: f64) -> (usize, usize, usize) {
        let mut ncell = [1usize; 3];

        for i in 0..3 {
            let len_i = self.sim_box.h.column(i).norm();
            ncell[i] = (len_i / rcut).floor() as usize;
            if ncell[i] == 0 {
                ncell[i] = 1;
            }
        }

        (ncell[0], ncell[1], ncell[2])
    }

    /// This function converts the numerical coordinate of cell into a single number.
    /// The number is the count of cells upto that cell in the simulation box
    pub fn cell_index(cx: usize, cy: usize, cz: usize, nx: usize, ny: usize) -> usize {
        (cz * ny + cy) * nx + cx
    }

    /// The parallel calculation of the lennard jones potential requires numerical cell coordinates
    /// rather than a single number. This function provides those coordinates.
    pub fn cell_indices_for_parallel(
        nx: usize,
        ny: usize,
        nz: usize,
    ) -> Vec<(usize, usize, usize)> {
        let mut indices: Vec<(usize, usize, usize)> = Vec::with_capacity(nx * ny * nz);

        for cx in 0..nx {
            for cy in 0..ny {
                for cz in 0..nz {
                    indices.push((cx, cy, cz));
                }
            }
        }
        indices
    }
}

/// This struct defines the coordinates of a cell
#[derive(Debug, Clone, Copy)]
pub struct Offset3 {
    pub dx: isize,
    pub dy: isize,
    pub dz: isize,
}

/// 14 forward neighbour offsets for a 3D cell-linked list:
pub const FORWARD_NEIGHBOUR_OFFSETS: [Offset3; 14] = [
    Offset3 {
        dx: 0,
        dy: 0,
        dz: 0,
    }, // self cell
    Offset3 {
        dx: 1,
        dy: 0,
        dz: 0,
    },
    Offset3 {
        dx: -1,
        dy: 1,
        dz: 0,
    },
    Offset3 {
        dx: 0,
        dy: 1,
        dz: 0,
    },
    Offset3 {
        dx: 1,
        dy: 1,
        dz: 0,
    },
    Offset3 {
        dx: -1,
        dy: -1,
        dz: 1,
    },
    Offset3 {
        dx: 0,
        dy: -1,
        dz: 1,
    },
    Offset3 {
        dx: 1,
        dy: -1,
        dz: 1,
    },
    Offset3 {
        dx: -1,
        dy: 0,
        dz: 1,
    },
    Offset3 {
        dx: 0,
        dy: 0,
        dz: 1,
    },
    Offset3 {
        dx: 1,
        dy: 0,
        dz: 1,
    },
    Offset3 {
        dx: -1,
        dy: 1,
        dz: 1,
    },
    Offset3 {
        dx: 0,
        dy: 1,
        dz: 1,
    },
    Offset3 {
        dx: 1,
        dy: 1,
        dz: 1,
    },
];
