#!/bin/bash

# If large system (more atoms/basis than in this example),
#   you may want to run these steps by submitting to a cluster

##### Run the TDDFT with additional keywords "IOp(6/8=3) IOp(9/40=5)" #####
#module load gaussian
g16 < geometry.com | tee geometry.out
formchk geometry.chk geometry.fchk

##### Run Multiwfn to get densities as cube files #####
../make_transition_density.slurm      | tee make_transition_density.out
../make_electric_dipole_density.slurm | tee make_electric_dipole_density.out
../make_magnetic_dipole_density.slurm | tee make_magnetic_dipole_density.out

##### Move all cube files to new folder #####
rm -r cube_files
mkdir cube_files
mv *cube cube_files/

# Run post-processing script to generate some plots
python3 ../compute_plots_projections.py