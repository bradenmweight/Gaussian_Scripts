#!/bin/bash

mkdir TDMat
cd TDMat
    cp ~/Gaussain_scripts/Transition_Density/make_rho_ij_low_density.slurm .
    sbatch make_rho_ij_low_density.slurm
    cd ../

mkdir NTOs
cd NTOs/
    cp ~/Gaussain_scripts/NTOs/submit.NTOs .
    sbatch submit.NTOs
    cd ../
