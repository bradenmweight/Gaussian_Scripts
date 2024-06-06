#!/bin/bash
#SBATCH -p debug
#SBATCH -J TransDens
#SBATCH -o output.slurm
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 1GB

cp ~/THUNDER_HOME/Multiwfn_3.6_bin_Linux_noGUI/Multiwfn_3.6_bin_Linux_noGUI/settings.ini .

# Generate ground state density: T_{00}
~/THUNDER_HOME/Multiwfn_3.6_bin_Linux_noGUI/Multiwfn_3.6_bin_Linux_noGUI/Multiwfn << EOF
../geometry.fchk
5
1
1
2
0
0
0
0
0
0
0
EOF

for state in {1..15}; do

echo "Generating first transition density for transition state ${state}."

~/THUNDER_HOME/Multiwfn_3.6_bin_Linux_noGUI/Multiwfn_3.6_bin_Linux_noGUI/Multiwfn << EOF
../geometry.fchk
18
1
../geometry.out
${state}
1
1
13
0
0
0
0
0
0
0
EOF

mv transdens.cub trans-0_${state}.cube

done
