#!/bin/bash

# Generate orbitals for the following G16 job:
###### %chk=geometry.chk
###### %mem=30GB
###### %nprocshared=12
######
###### #P wB97XD/6-31G guess=read POP=NBORead
######
###### Title Card Required
######
###### 0 1
######
###### $NBO plot $END




PATH_TO_MULTIWFN="/home/bweight/Multiwfn_3.7_bin_Linux_noGUI/Multiwfn"
HOMO_LABEL=$(grep 'alpha electrons' GS.out | awk '{print $1}') # For H-F system (BLYP/6-31G*), the HOMO orbital was #5 since there were 10 electrons
GRID_QUALITY=2 # 1 = Coarse, 2 = Medium, 3 = Fine # For H-F system, "2" quality produced 7.3 MB files

# OBTAIN "HOMO_NBO" CUBE FILE
${PATH_TO_MULTIWFN} << EOF
FILE.31
FILE.37
5
4
${HOMO_LABEL}
${GRID_QUALITY}
2
0
0
0
0
0
0
EOF

mv MOvalue.cub HOMO_NBO.cube


# OBTAIN "LUMO_NBO" CUBE FILE
${PATH_TO_MULTIWFN} << EOF
FILE.31
FILE.37
5
4
$((HOMO_LABEL+1))
${GRID_QUALITY}
2
0
0
0
0
0
0
EOF

mv MOvalue.cub LUMO_NBO.cube
