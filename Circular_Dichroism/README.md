# Circular Dichroism
    ~~~ Braden M. Weight -- May 9, 2024 ~~~

## Theory
Electric Dipole Moment:
$\hat{\mu}^\mathrm{el} = -\sum_{j}^{N_\mathrm{el}} {\bf \hat{r}}_j + \sum_{J}^{N_\mathrm{ions}} Z_J {\bf R}_J \otimes \mathbb{I}_\mathrm{el}$

Magnetic Dipole Moment:
$\hat{\mu}^\mathrm{mag} = -\sum_{j}^{N_\mathrm{el}} {\bf \hat{r}}_j \times \frac{i}{\hbar}{\bf \hat{p}}_j $

Note: $\hat{p}_x = -i\hbar \hat{\nabla}_x$

Circular Dichroism Oscillator Strength: $$f_{0\alpha}^\mathrm{CD} = 2 E_{0\alpha}~\hat{\mu}^\mathrm{el} \cdot \hat{\mu}^\mathrm{mag},$$ where $E_{0\alpha} = E_\alpha - E_0$

Circular Dichroism Spectroscopy: $$\mathrm{CD}(E) = \sum_\alpha f_{0\alpha}^\mathrm{CD} \delta(E - E_{0\alpha}) \approx \sum_\alpha f_{0\alpha}^\mathrm{CD} \mathrm{e}^{-\frac{(E - E_{0\alpha})^2}{2\sigma^2}},$$ where $\sigma$ is a finite broadening parameter

Above, 
$$\boldsymbol{\mu}_{0\alpha}^\mathrm{el} = \mathrm{Tr}[\hat{\mu}^\mathrm{el}~\hat{\xi}_{0\alpha}] = \int dr \mu^\mathrm{el}(r)~\xi_{0\alpha}(r) = \int d\boldsymbol{r}~\psi^*_\alpha(\boldsymbol{r})~\boldsymbol{r}~\psi_0(\boldsymbol{r}')$$ 
and 
$$\mu_{0\alpha}^\mathrm{mag} = \mathrm{Tr}[\hat{\mu}^\mathrm{mag}~\hat{\xi}_{0\alpha}] = \int d\boldsymbol{r} \mu^\mathrm{mag}(\boldsymbol{r})~\xi_{0\alpha}(\boldsymbol{r}) = \int d\boldsymbol{r}~\psi^*_\alpha(\boldsymbol{r})~(\boldsymbol{r}\times \boldsymbol{\hat{\nabla}})~\psi_0(\boldsymbol{r})$$
where $\xi_{0\alpha}(r,r') = \psi^*_\alpha(r)\psi_0(r')$ is the one-particle transition density matrix between the ground $0$ and excited state $\alpha$.


## Intructions
```
1. Gaussian TDDFT with "IOp(6/8=3) IOp(9/40=5)"
2. make_dipole_density.slurm
3. make_magnetic_dipole_density.slurm
4. compute_CD_Density.py
```
Note: This script computes the rotary strength density from the transition dipole and transition magnetic dipole densities as output from MultiWfn (V3.7).