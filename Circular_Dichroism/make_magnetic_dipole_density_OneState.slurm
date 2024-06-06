#!/bin/bash
#SBATCH -p exciton -A exciton
#SBATCH -J TransDens
#SBATCH -o output.slurm
#SBATCH -t 5-00:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --mem 30GB

state=10 # Number of excited states
GRID=3 # 1 for coarse grid, 2 for fine grid, 3 for ultrafine grid

cp ~/Multiwfn_3.7_bin_Linux_noGUI/settings.ini .

# Generate X-projected transition dipole density
~/Multiwfn_3.7_bin_Linux_noGUI/Multiwfn << EOF
../geometry.fchk
18
1
../geometry.out
${state}
-1
1
${GRID}
17
1
0
0
0
0
0
0
0
EOF

mv magtrdipdens.cub MAG_DIP_DENSITY_0_${state}_X.cube


# Generate Y-projected transition dipole density
~/Multiwfn_3.7_bin_Linux_noGUI/Multiwfn << EOF
../geometry.fchk
18
1
../geometry.out
${state}
-1
1
${GRID}
17
2
0
0
0
0
0
0
0
EOF

mv magtrdipdens.cub MAG_DIP_DENSITY_0_${state}_Y.cube


# Generate Z-projected transition dipole density
~/Multiwfn_3.7_bin_Linux_noGUI/Multiwfn << EOF
../geometry.fchk
18
1
../geometry.out
${state}
-1
1
${GRID}
17
3
0
0
0
0
0
0
0
EOF

mv magtrdipdens.cub MAG_DIP_DENSITY_0_${state}_Z.cube
