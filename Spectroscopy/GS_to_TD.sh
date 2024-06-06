#####!/bin/bash

mydir=$(pwd)
for d in $(find . -type d -name "TD-singlets"); do

    if [[ $d == *'7'* ]]; then
        echo "I Skipped #7: $d/"
        continue
    fi
    echo "$d"

    cd $d/;
    #python ~/Gaussain_scripts/Excited_States/GS_to_TD.py
    cp ~/submit_scripts/submit.gaussian .
    sbatch submit.gaussian
    cd $mydir/

done




