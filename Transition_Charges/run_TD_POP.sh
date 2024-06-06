#!/bin/bash

# Run this code inside a TD_DFT/TD_CHARGES directory
# Need GS.out, GS.fchk, TD.out, and TD.fchk
# This will generate transition charges Q_{0k} and Q_{kk}
# I don't know how to get Q_{j!=k} yet
###### February 3, 2023 -- Braden M. Weight ######



TEMPLATE_FILE="geometry.com_TEMPLATE" # Template TD file with root=BASH_ROOT in TD command
                                      # TD=(read,singlets,nstates=100,root=BASH_ROOT) pop=Full
PATH_TO_MULTIWFN="/home/bweight/Multiwfn_3.7_bin_Linux_noGUI/Multiwfn"


### Use Gaussian to get diagonal charges ###

rm -r TD
mkdir TD

for state in {1..100}; do
    cd TD/
        mkdir S${state}
        cd S${state}/
            cp ../../${TEMPLATE_FILE} geometry.com
            sed -i "s/BASH_ROOT/$state/g" geometry.com
            cp ../../submit.gaussian .
            #sbatch submit.gaussian
            cd ../../
done

### Use MultiWfn to get g-to-e transition charges ###

rm -r DATA
mkdir DATA

# OBTAIN GS Atomic Charges
${PATH_TO_MULTIWFN} << EOF
GS.fchk
7
5
1
y
EOF

cat atmtrchg.chg | awk '{print $5}' > DATA/Q_0-0.dat

for state in {1..100}; do
# OBTAIN Atomic Transition Charges
${PATH_TO_MULTIWFN} << EOF
TD.fchk
18
12
TD.out
${state}
EOF

cat atmtrchg.chg | awk '{print $5}' > DATA/Q_0-${state}.dat

done



### We are missing e-to-e transition charges ###

