# Circular Dichroism
    ~~~ Braden M. Weight -- May 9, 2024 ~~~

## Theory
```
\hat{\mu}^\mathrm{el} = -\sum_{i}^{N_\mathrm{el}} \hat{\bf r} + \sum_{I}^{N_\mathrm{ions}} \hat{\bf R}
```



    Intructions:
    1. Gaussian TDDFT with Additional Keywords: "IOp(6/8=3) IOp(9/40=5)"
    2. make_dipole_density.slurm
    3. make_magnetic_dipole_density.slurm
    4. compute_CD_Density.py

Note: This script computes the rotary strength density from the transition dipole and transition magnetic dipole densities as output from MultiWfn (V3.7).