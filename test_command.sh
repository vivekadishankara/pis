# This script will run with the example lammps data file in the test folder with the starting temperature as 5K for 100 steps
time RAYON_NUM_THREADS=4 cargo run --release --bin pis --features cli -- -i example/input.pis
