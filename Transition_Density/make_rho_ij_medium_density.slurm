#!/bin/bash
#SBATCH -p debug
#SBATCH -J TransDens
#SBATCH -o output.slurm
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 1GB

cd $PBS_O_WORKDIR

cp ~/Multiwfn_3.5_bin_Linux_noGUI_GIFF/settings.ini .

for n in {1..15}; do

~/Multiwfn_3.5_bin_Linux_noGUI_GIFF/Multiwfn << EOF
geometry.fchk
18
1
geometry.out
${n}
1
2
13
0
0
0
0
0
0
0
EOF

mv transdens.cub trans-${n}.cube

done

















