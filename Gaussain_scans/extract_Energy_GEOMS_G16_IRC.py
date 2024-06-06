import numpy as np
from matplotlib import pyplot as plt
import subprocess as sp
import os

DO_SUBMIT_TDDFT = True
DO_QCHEM        = True # QCHEM or G16

def extract_IRC_SCAN():
    global NGEOMS, NATOMS
    NGEOMS        = int(sp.check_output( "grep 'Point Number:' geometry.out | wc -l", shell=True ) )
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
        #if ( "Optimization completed" in line ): # Do not put a period here.
        if ( len(line.split()) > 2 ):
            if ( line.split()[:2] == "Point Number:".split() or line.split(":")[0].strip() == "Point Number:".split(":")[0].strip() ):
                print( count, int(line.split(":")[1].split()[0]), len(ENERGY_SCAN) )
                if ( int(line.split(":")[1].split()[0]) == len(ENERGY_SCAN)+1 ):
                    ENERGY_SCAN.append( ENERGY_TMP )
                    GEOMS.append( GEOM_TMP )
    # Catch final geometry and energy
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

    return LABELS, GEOMS

def make_input_qchem( LABELS, GEOM ):
    string = """\
$molecule
   0 1
"""

    for at in range( NATOMS ):
        string += f"   {LABELS[at]} {' '.join(GEOM[at].astype(str))}\n"
    string += "$end\n"
    
    string += """\
   $rem
   METHOD             wB97XD
   RPA                TRUE
   BASIS              6-31G*
   CIS_N_ROOTS        30
   CIS_SINGLETS       TRUE
   CIS_TRIPLETS       FALSE
   MAKE_CUBE_FILES    TRUE
   PLOTS              TRUE
   STS_MOM            TRUE
   SYM_IGNORE         TRUE
   SYMMETRY           FALSE
   MEM_TOTAL          2000
$end

$plots
   grid_points                    75 75 75
   total_density                  0-30
   transition_density             1-30
$end
"""

    FILE01 = open("QCHEM.in","w")
    FILE01.write( string )
    FILE01.close()


def main():
    LABELS, GEOMS = extract_IRC_SCAN()

    if ( DO_SUBMIT_TDDFT ):
        sp.call("rm -r TDDFT", shell=True)
        sp.call("mkdir TDDFT", shell=True)
        os.chdir("TDDFT/")
        for step in range( NGEOMS ):
            print(f"TDDFT for geometry {step}")
            sp.call(f"mkdir STEP_{step}", shell=True)
            os.chdir(f"STEP_{step}/")
            if ( DO_QCHEM ):
                make_input_qchem( LABELS, GEOMS[step] )
                sp.call("cp ~/submit_scripts/submit.QCHEM .", shell=True)
                sp.call("sbatch submit.QCHEM", shell=True)

            else: # TODO
                print("G16 not yet implemented. Use QCHEM.")
                # make_input_g16( LABELS, GEOMS )
            os.chdir("../")
        os.chdir("../")


if ( __name__ == "__main__" ):
    main()