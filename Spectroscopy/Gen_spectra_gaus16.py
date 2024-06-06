# SYNTAX: python3 Gen_spectra_gaus16.py

import numpy as np
import subprocess as sp
from matplotlib import pyplot as plt

# For plotting
PLOT_DOS  = False
do_CD_LEN = False # In addition to velocity gauge, show length gauge
NPTS = 1000
SIG  = 0.025   # eV
EMIN = 1.0   # eV
EMAX = 2.0   # eV

###########################################################

E_GRID = np.linspace( EMIN, EMAX, NPTS )
dE     = E_GRID[1] - E_GRID[0]

# Get data from Gaussian
sp.call("rm EXC.dat OSC.dat CHECK_EOMCCSD.dat",shell=True)
sp.call("grep 'EOMCCSD' geometry.out -i | awk '{print $5}' > CHECK_EOMCCSD.dat",shell=True)
sp.call("grep 'Excited State' geometry.out | awk '{print $5}' > EXC.dat",shell=True)
sp.call("grep 'Excited State' geometry.out | awk '{print $9}' > OSC.dat",shell=True)
Energies = np.loadtxt("EXC.dat")
print(Energies)
OSC_STR  = np.array([j[2:] for j in open("OSC.dat","r").readlines()]).astype(float)
if ( len(open("CHECK_EOMCCSD.dat","r").readlines()) > 0 ):
  print("Found EOM-CCSD")
  NCIS = int(sp.check_output("grep 'Convergence achieved on expansion vectors' geometry.out -B 1 | head -n 1 | awk '{print $2}'",shell=True))
  #NCIS = int(sp.check_output("grep 'initial guesses have been' geometry.out | awk '{print $1}'",shell=True))
  NSTATES  = len(Energies) - NCIS
  Energies = Energies[NCIS:]
  OSC_STR  = OSC_STR[NCIS:]
else:
  print("Found TD-DFT")
  NSTATES  = len(Energies)
sp.call("rm EXC.dat OSC.dat CHECK_EOMCCSD.dat",shell=True)
np.savetxt("Energy_Osc.dat", np.c_[Energies, OSC_STR], fmt="%1.6f")

# Get Rotary Strengths
sp.call("rm ROT_VEL.dat ROT_LEN.dat",shell=True)
sp.call(f"grep 'R(velocity)' geometry.out -A {NSTATES} | tail -n {NSTATES} "+" | awk '{print $5}' > ROT_VEL.dat",shell=True)
sp.call(f"grep 'R(length)' geometry.out -A {NSTATES} | tail -n {NSTATES} "+"| awk '{print $5}' > ROT_LEN.dat",shell=True)
sp.call(f"grep 'transition magnetic dipole moments (Au):' geometry.out -A {NSTATES+1} | tail -n {NSTATES} "+"| awk '{print $2, $3, $4}' > MAG_DIP.dat",shell=True)
sp.call(f"grep 'transition electric dipole moments (Au):' geometry.out -A {NSTATES+1} | tail -n {NSTATES} "+"| awk '{print $2, $3, $4}' > EL_DIP.dat",shell=True)
MAG_DIP = np.loadtxt("MAG_DIP.dat")
EL_DIP = np.loadtxt("EL_DIP.dat")
ROT_VEL = np.loadtxt("ROT_VEL.dat")
ROT_LEN = np.loadtxt("ROT_LEN.dat")
sp.call("rm ROT_VEL.dat ROT_LEN.dat EL_DIP.dat MAG_DIP.dat",shell=True)
np.savetxt("ROTARY_STRENGTHS.dat", np.c_[Energies, ROT_VEL, ROT_LEN], fmt="%1.6f", header="E(ev) Rotary(vel) Rotary(len)")
np.savetxt("MAGNETIC_DIPOLE.dat", np.c_[Energies, MAG_DIP[:,0], MAG_DIP[:,1], MAG_DIP[:,2]], fmt="%1.6f", header="E(ev) Mag.Dipole(a.u.)X,Y,Z")
np.savetxt("ELECTRIC_DIPOLE.dat", np.c_[Energies, EL_DIP[:,0], EL_DIP[:,1], EL_DIP[:,2]], fmt="%1.6f", header="E(ev) El.Dipole(a.u.)X,Y,Z")

SPEC   = np.zeros( (NPTS) )
CD_VEL = np.zeros( (NPTS) )
CD_LEN = np.zeros( (NPTS) )
DOS    = np.zeros( (NPTS) )
for pt in range(NPTS):
  E = E_GRID[pt]
  SPEC[pt]   = np.sum( OSC_STR[:] * np.exp( -(E-Energies[:])**2 / 2 / SIG**2) )
  DOS[pt]    = np.sum( 1.00000000 * np.exp( -(E-Energies[:])**2 / 2 / SIG**2) )
  CD_VEL[pt] = np.sum( ROT_VEL[:] * np.exp( -(E-Energies[:])**2 / 2 / SIG**2) )
  CD_LEN[pt] = np.sum( ROT_LEN[:] * np.exp( -(E-Energies[:])**2 / 2 / SIG**2) )

np.savetxt(f'SPEC__SIG_{SIG}_NStates_{NSTATES}.dat', np.c_[E_GRID, SPEC], fmt="%1.6f")
np.savetxt(f'DOS__SIG_{SIG}_NStates_{NSTATES}.dat',  np.c_[E_GRID, DOS],  fmt="%1.6f")
np.savetxt(f'CD_VEL__SIG_{SIG}_NStates_{NSTATES}.dat',  np.c_[E_GRID, CD_VEL],  fmt="%1.6f")
np.savetxt(f'CD_LEN__SIG_{SIG}_NStates_{NSTATES}.dat',  np.c_[E_GRID, CD_LEN],  fmt="%1.6f")

# Simple Plot
plt.plot( E_GRID, SPEC, c='blue', label="Abs." )
if ( PLOT_DOS == True ):
  plt.plot( E_GRID, DOS / np.max(DOS) * np.max(SPEC), c='blue',alpha=0.25, linewidth=6, label="DOS" )
#plt.scatter( Energies, OSC_STR, c='blue' )
markerline, stemlines, baseline = plt.stem(Energies, OSC_STR, linefmt="blue", markerfmt="o", )
markerline.set_markerfacecolor('none')
markerline.set_markeredgecolor('blue')
markerline.set_markersize(8)
markerline.set_markeredgewidth(1.5)
plt.xlabel('Energy (eV)', fontsize=15)
plt.ylabel('Absorption (Osc. Str.)', fontsize=15)
plt.xlim(EMIN,EMAX)
plt.ylim(0)
plt.legend()
plt.tight_layout()
plt.savefig(f"SPEC__SIG_{SIG}_NStates_{NSTATES}.jpg", dpi=300)
plt.clf()

# Simple Plot
plt.plot( E_GRID, CD_VEL, "-",  c='blue', label="CD (Vel.)" )
if ( do_CD_LEN == True ):
  plt.plot( E_GRID, CD_LEN, "-", c='red', label="CD (Len.)" )
if ( PLOT_DOS == True ):
  plt.plot( E_GRID, DOS, c='blue',alpha=0.25, linewidth=6, label="DOS" )
#plt.scatter( Energies, OSC_STR, c='blue' )
markerline, stemlines, baseline = plt.stem(Energies, ROT_VEL, linefmt="red", markerfmt="o", )
markerline.set_markerfacecolor('none')
markerline.set_markeredgecolor('red')
markerline.set_markersize(8)
markerline.set_markeredgewidth(1.5)
if ( do_CD_LEN == True ):
  markerline, stemlines, baseline = plt.stem(Energies, ROT_LEN, linefmt="blue", markerfmt="o", )
  markerline.set_markerfacecolor('none')
  markerline.set_markeredgecolor('blue')
  markerline.set_markersize(8)
  markerline.set_markeredgewidth(1.5)
plt.xlabel('Energy (eV)', fontsize=15)
plt.ylabel('CD (10**-40 erg-esu-cm/Gauss)', fontsize=15)
plt.xlim(EMIN,EMAX)
#plt.ylim(0)
plt.legend()
plt.tight_layout()
plt.savefig(f"CD__SIG_{SIG}_NStates_{NSTATES}.jpg", dpi=300)
plt.clf()











