#!/bin/bash
#SBATCH -p debug
#SBATCH -J TransDens
#SBATCH -o output.slurm
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --mem 30GB

GRID=2 # 1 for coarse grid, 2 for fine grid, 3 for ultrafine grid
OUTFILE="geometry.out"
FCHKFILE="geometry.fchk"
MULTIWFN="/home/bweight/Multiwfn_3.7_bin_Linux_noGUI/Multiwfn"

#################################################################

NSTATES=$(grep "Excited State" ${OUTFILE} | tail -n 1 | awk '{print $3}')
NSTATES=${NSTATES:0:1}
echo "!!!!   I found ${NSTATES} excited state(s) in ${OUTFILE}."

# Generate transition density
for ((state = 1 ; state < $NSTATES + 1 ; state++ )); do

# If only one state in the calculation, set state to blank
# This is because Multiwfn will skip this input line if only one state...
# Else set the state to the current state
statemultiwfn=${state}
if [ $NSTATES -eq 1 ]; then
    statemultiwfn=""
fi

echo "Generating first transition density for transition state ${state}."

${MULTIWFN} << EOF
${FCHKFILE}
18
1
${OUTFILE}
${statemultiwfn}
1
${GRID}
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
