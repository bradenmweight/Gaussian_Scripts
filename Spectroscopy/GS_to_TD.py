# Read .out (.log) file from AM1 calc. and obtain newest geoemtry during optimization
#   make CAM-B3LYP jobs starting from most recent geometry.

# SYNTAX: python3 AM1_output_to_CAM-B3LYP_input.py

import subprocess as sp
import glob, sys, os
import numpy as np

# Get geometry
if ( os.path.isfile("geometry.out") ):
        name = "geometry.out"
elif ( os.path.isfile("geometry.out-old" ) ):
        name = "geometry.out-old"
else:
        print ("Output file not found. Quitting.")
        exit()
sp.call("module load openbabel", shell=True)
sp.call(f"/software/openbabel/2.4.1/bin/babel -iout {name} -oxyz geometry.xyz", shell=True)
XYZLines = open("geometry.xyz","r").readlines()
Coords = np.array ( [ XYZLines[j].split() for j in range(len(XYZLines)) if len(XYZLines[j].split()) == 4 ] )
unique_atoms = []
for n in Coords:
	if ( n[0] not in unique_atoms ):
		unique_atoms.append(n[0])

print ("Unique Atom Types", unique_atoms)

# Clean files
sp.call("rm Gau-* e *.e*",shell=True)
sp.call("mv geometry.out geometry.out-old",shell=True)

# Copy PBS file to CAM directory
sp.call("cp ~/submit_scripts/submit.gaussian .", shell=True)

# Make new geometry.com file only using geometry from AM1 calc.
newComFile = open("geometry.com","w")

newComFile.write("%chk=geometry.chk\n")
newComFile.write("%mem=139GB\n")
newComFile.write("%nprocshared=8\n\n")

#if ("Cl" in unique_atoms or "O" in unique_atoms or "N" in unique_atoms):
#        newComFile.write("#p SCF=XQC OPT CAM-B3LYP/GEN\n\n")
#else:
newComFile.write("#p TD=(singlets,nstates=15) cam-b3lyp/sto-3g\n#p IOp(6/8=3) IOp(9/40=4)\n\n")

newComFile.write("TitleMe\n\n")

newComFile.write("0 1\n")

for n in range(len(Coords)):
	outArray = [ Coords[n,0], Coords[n,1], Coords[n,2], Coords[n,3]  ]
	newComFile.write(" ".join(outArray) + "\n")

# IF necessary, write mixed basis
#if ( "Cl" not in unique_atoms and "O" not in unique_atoms and "N" not in unique_atoms):
newComFile.write("\n\n\n\n\n\n\n\n\n\n\n")
"""
else:
        unique_atoms.remove("C")
        unique_atoms.remove("H")
        newComFile.write("\n")
        newComFile.write(" ".join(unique_atoms) + " 0\n")
        newComFile.write("6-31G*\n")
        newComFile.write("****\n")
        newComFile.write("C H 0\n")
        newComFile.write("STO-3G\n")
        newComFile.write("****\n\n\n\n\n\n\n\n\n\n")
"""

sp.call("sbatch submit.gaussian",shell=True)

os.chdir("../")


