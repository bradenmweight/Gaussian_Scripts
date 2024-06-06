import numpy as np
from matplotlib import pyplot as plt
import subprocess as sp

NGEOMS        = int(sp.check_output( "grep 'on scan point' geometry.out | tail -n 1 | awk '{print $16}'", shell=True ) )
NATOMS        = int(sp.check_output( "grep 'NAtoms=' geometry.out | head -n 1 | awk '{print $2}'", shell=True ) )
LABELS, GEOMS = [], []
ENERGY_SCAN   = []

print( "Number of Atoms: ", NATOMS )
print( "Number of Geometries: ", NGEOMS )

FILE01 = open("geometry.out","r").readlines()
for count,line in enumerate( FILE01 ):
    if ( "SCF Done" in line ):
        ENERGY_TMP = float( line.split()[4] )
    if ( "Input orientation:" in line ):
        GEOM_TMP  = np.zeros( (NATOMS,3) )
        lines = FILE01[count+5:count+5+NATOMS]
        for at in range( NATOMS ):
            LABELS.append( lines[at].split()[1] )
            GEOM_TMP[at,:] = np.array( lines[at].split()[3:] ).astype(float)
    if ( "Optimization completed" in line ): # Do not put a period here.
        ENERGY_SCAN.append( ENERGY_TMP )
        GEOMS.append( GEOM_TMP )

GEOMS        = np.array( GEOMS )
ENERGY_SCAN  = np.array( ENERGY_SCAN ) * 627.5
LABELS       = np.array(LABELS[:NATOMS]).astype(int)

plt.plot( np.arange(len(ENERGY_SCAN)), ENERGY_SCAN - ENERGY_SCAN[0], "-o", c='black' )
plt.xlabel( "Scan Step", fontsize=15 )
plt.ylabel( "Energy (kcal/mol)", fontsize=15 )
plt.savefig( "energy_scan.jpg", dpi=300 )
plt.clf()
np.savetxt( "energy_scan.dat", np.c_[np.arange(len(ENERGY_SCAN)), ENERGY_SCAN], fmt="%1.0f %1.4f", header="Step   Energy(kcal/mol)")
np.savetxt( "energy_scan_reversed.dat", np.c_[np.arange(len(ENERGY_SCAN)), ENERGY_SCAN[::-1]], fmt="%1.0f %1.4f", header="Step   Energy(kcal/mol)")
plt.clf()


# Make XYZ file
DICT = \
{    1: "H",   2: "He",  3: "Li",  4: "Be",  5: "B",   6: "C",   7: "N",   8: "O",   9: "F",  10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",  16: "S",  17: "Cl", 18: "Ar", 19: "K",  20: "Ca",
    21: "Sc", 22: "Ti", 23: "V",  24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn",
    31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br"}
FILE01 = open("trajectory.xyz","w")
for geom in range( NGEOMS ):
    FILE01.write( str(NATOMS) + "\n" )
    FILE01.write( "Scan Step " + str(geom) + "\n" )
    try:
        for at in range( NATOMS ):
            FILE01.write( DICT[LABELS[at]] + " " + " ".join( GEOMS[geom,at].astype(str) ) + "\n" )
    except IndexError:
        print("Warning!!! Probably the scan did not complete. Be careful.")
        break

FILE01.close()

FILE01 = open("trajectory_reversed.xyz","w")
for geom in range( NGEOMS-1, -1, -1):
    FILE01.write( str(NATOMS) + "\n" )
    FILE01.write( "Scan Step " + str(geom) + "\n" )
    try:
        for at in range( NATOMS ):
            FILE01.write( DICT[LABELS[at]] + " " + " ".join( GEOMS[geom,at].astype(str) ) + "\n" )
    except IndexError:
        print("Warning!!! Probably the scan did not complete. Be careful.")
        break
FILE01.close()