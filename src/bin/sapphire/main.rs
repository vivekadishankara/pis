use std::{fs::File, io::{BufWriter, Write}};

use nalgebra::{Matrix3, Matrix3xX, Vector3};

fn main() {
    generate_fcc_supercell(5.41, 5, 5, 5, "example.txt").expect("Generating the file led to an error")
}

fn fcc_h_matrix(a: f64) -> Matrix3<f64> {
    Matrix3::identity() * a
}

fn fcc_basis_frac() -> Matrix3xX<f64> {
    let a1: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
    let a2: Vector3<f64> = Vector3::new(0.0, 0.5, 0.5);
    let a3: Vector3<f64> = Vector3::new(0.5, 0.0, 0.5);
    let a4: Vector3<f64> = Vector3::new(0.5, 0.5, 0.0);

    Matrix3xX::from_columns(&[a1, a2,a3, a4])
}

fn generate_super_cell(
    h: Matrix3<f64>, 
    basis_frac: Matrix3xX<f64>, 
    nx: usize, ny: usize, nz: usize,
    path: &str,
) -> std::io::Result<()> {
    let n_atoms = nx * ny * nz * basis_frac.ncols();
    let mut ncols: usize = 0; 

    let file = File::create(path)?;
    let mut out  = BufWriter::new(file);

    if path.contains(".xyz") {
        writeln!(out, "{}", n_atoms)?;
        writeln!(out, "Generated cell")?;
    }
    

    if path.contains(".txt") {
        writeln!(out, "{} atoms", n_atoms)?;
        writeln!(out, "1 atom types\n")?;

        writeln!(out, "0.0 {} xlo xhi", nx as f64 * h[(0,0)])?;
        writeln!(out, "0.0 {} ylo yhi", ny as f64 * h[(0,0)])?;
        writeln!(out, "0.0 {} zlo zhi\n", nz as f64 * h[(0,0)])?;

        writeln!(out, "Masses\n1 39.948\n")?;
        writeln!(out, "PairCoeffs\n1 0.23 3.405 8.5\n")?;

        writeln!(out, "Atoms")?;
    }

    let a_type = "I";

    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let cell_origin_frac: Vector3<f64> = Vector3::new(ix as f64, iy as f64, iz as f64);

                for b in basis_frac.column_iter() {
                    let frac = cell_origin_frac + b;

                    let pos = &h * frac;
                    ncols += 1;
                    if path.contains(".xyz") {
                        writeln!(out, "{} {} {} {}", a_type, pos[0], pos[1], pos[2])?;
                    } else {
                        writeln!(out, "{} 1 {} {} {}", ncols, pos[0], pos[1], pos[2])?;
                    }
                    
                }
            }
        }
    }
    Ok(())
}

fn generate_fcc_supercell(a: f64, nx: usize, ny: usize, nz: usize, path: &str) -> std::io::Result<()> {
    let h = fcc_h_matrix(a);
    let basis_frac = fcc_basis_frac();
    generate_super_cell(h, basis_frac, nx, ny, nz, path)?;

    Ok(())
}
