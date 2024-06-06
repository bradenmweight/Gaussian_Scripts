This script computes the rotary strength density from the 
    transition dipole and transition magnetic dipole densities
    as output from MultiWfn (V3.7).
    
    ~~~ Braden M. Weight -- May 9, 2024 ~~~

    Intructions:
    1. Gaussian TDDFT with Additional Keywords: "IOp(6/8=3) IOp(9/40=5)"
    2. make_dipole_density.slurm
    3. make_magnetic_dipole_density.slurm
    4. compute_CD_Density.py
