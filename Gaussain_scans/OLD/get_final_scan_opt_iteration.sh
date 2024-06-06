#!/bin/bash

N=26 # Max number of iterations = Scan Steps + 1

# May need to grep "out of" to get spacing right. Sometimes, it's weird. Then replace
#  first number with ${n} and last number with ${N}.
# Also, change length of for loop to be ${N. Couldn't figure out how to hard code this...

for n in {1..26}; do line=$(grep " ${n} out of    ${N}" geometry.out | tail -n 1); echo $line; echo $line >> Scan_steps_Final.dat; done
#for n in {1..101}; do (grep " ${n} out of   ${N}" geometry.out | tail -n 1) >> Scan_steps_Final.dat; done
