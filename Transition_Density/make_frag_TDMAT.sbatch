#!/bin/bash
#SBATCH -p preempt
#SBATCH -J TransDens
#SBATCH -o output.slurm
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem 5GB

module load gaussian
module load openbabel

/software/openbabel/2.4.1/bin/babel -iout ../geometry.out -ocom geometry.com_FULL

cp ~/THUNDER_HOME/transition_density_scripts/slice_tube_output_fragments.py .
python3 slice_tube_output_fragments.py

cp ~/THUNDER_HOME/Multiwfn_3.6_bin_Linux_noGUI/Multiwfn_3.6_bin_Linux_noGUI/settings.ini .
for n in {1..15}; do

echo "Generating atomic TDM for transition state 1."
~/THUNDER_HOME/Multiwfn_3.6_bin_Linux_noGUI/Multiwfn_3.6_bin_Linux_noGUI/Multiwfn << EOF
../geometry.fchk
18
9
1
../geometry.out
${n}
0
n
EOF

~/THUNDER_HOME/Multiwfn_3.6_bin_Linux_noGUI/Multiwfn_3.6_bin_Linux_noGUI/Multiwfn << EOF
../geometry.fchk
18
2
tdmat.txt
2
-1
0
fragdef.txt
3
EOF

#mv tdmat.txt tdmat_G_S${n}.txt
rm tdmat.txt # Removing this save lots of disk space
mv tmat.txt tmat_G_S${n}.txt

python3 ~/THUNDER_HOME/transition_density_scripts/plot_CNT_slice_TDMAT.py tmat_G_S${n}.txt
python3 ~/THUNDER_HOME/transition_density_scripts/calc_Ld_Lc_from_TDMAT.py tmat_G_S${n}.txt

done


