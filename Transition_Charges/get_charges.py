import numpy as np
import subprocess as sp
import os

# NEED TO RUN "run_TD_POP.sh" code



NSTATES = 100 # Excited States
CHARGES = []

for state in range( NSTATES ):
    print(state+1)
    os.chdir(f"TD/S{state+1}/")
    if ( state == 0 ):
        NATOMS = int( sp.check_output("grep 'NAtoms=' geometry.out | head -n 1 | awk '{print $2}'",shell=True) )
        CHARGES = np.zeros(( NSTATES, NATOMS ))
    sp.call("grep 'Mulliken charges:' geometry.out -A "+f"{NATOMS+1}"+" | tail -n "+f" {NATOMS}"+" | awk '{print $3}' > CHARGES.dat",shell=True)
    print(np.loadtxt("CHARGES.dat"))
    CHARGES[state,:] = np.loadtxt("CHARGES.dat")

    os.chdir("../../")

np.save("CHARGES.dat.npy", CHARGES)


