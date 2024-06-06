#!/bin/bash
#SBATCH -p exciton -A exciton
#SBATCH -J TransDens
#SBATCH -o output.slurm
#SBATCH -t 5-00:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --mem 30GB

NSTATES=10 # Number of excited states

cp ~/Multiwfn_3.7_bin_Linux_noGUI/settings.ini .

# Generate transition density
for ((state = 0 ; state < $NSTATES + 1 ; state++ )); do

echo "Generating first transition density for transition state ${state}."

~/Multiwfn_3.7_bin_Linux_noGUI/Multiwfn << EOF
../geometry.fchk
18
1
../geometry.out
${state}
1
3
13
0
0
0
0
0
0
0
EOF

mv transdens.cub TRANSITION_DENSITY_0_${state}.cube

done
