# The purpose of this script is to calculate Ld and Lc
#  from the transition density matrix, as described in 
#  Chem. Rev. 2002, 102, 3171-3212 by Tretiak/Mukamel.

# SYNTAX: python3 calc_Ld_Lc_from_TDMAT.py tmat_G_S1.txt
# (output from multiwfn 3.6 options: 18-->2)

import numpy as np
import sys

file05 = open("slice_width.txt","r")
dL = float(file05.readlines()[0])
file05.close()
print ("Slice Width:", dL)


name = sys.argv[1]
file01 = open(name,'r')
lines = file01.readlines()
N = int(len(lines) ** (1/2))
print ("Dimension N of Matrix:", N)

TDMAT = np.zeros((N,N))

file01.close()
file01 = open(name,'r')

for count, line in enumerate(file01):
  s = line.split()
  if (len(s) > 0):
    n = int(s[0])-1
    m = int(s[1])-1
    TD = float(s[2])
    TDMAT[n,m] = TD

# CALCULATE NORMS FOR PROBABILITIES
SUM1 = 0 # Sum of diagonal elements, Pn_DENOM
SUM2 = 0 # Sum of all elements, Qnm_DENOM
for i in range(N):
  for j in range(N):
    if (i == j):
      SUM1 += np.abs(TDMAT[i,i])
    SUM2 += np.abs(TDMAT[i,j])

SUM12 = SUM1*SUM1 # Square of diag sum
SUM22 = SUM2*SUM2 # Square of all sum

# CALCULATE Ld ~~ Size of Exciton (Localization)
Ld_DENOM = 0
for i in range(N):
  Ld_DENOM += TDMAT[i,i] * TDMAT[i,i]

Ld = SUM12 / Ld_DENOM


# CALCULATE Lc ~~ e-h Binding Length (Coherence)
Lc_DENOM = 0
for i in range(N):
  for j in range(N):
    Lc_DENOM += TDMAT[i,j] * TDMAT[i,j]

Lc = (1/Ld) * SUM22 / Lc_DENOM

Ld = Ld * dL / 10 # Convert to nm from index
Lc = Lc * dL / 10 # Convert to nm from index
print ("\tLd =", Ld, "nm \n\t\tMultiplied by dL for length.")
print ("\tLc =", Lc, "nm \n\t\tMultiplied by dL for length.")
print ("Check Lc =", (dL/10) * SUM22/SUM12 * Ld_DENOM / Lc_DENOM)

file02 = open(name + "_Ld_Lc.dat","w")
file02.write("Ld\tLc\n")
file02.write(str(Ld) + "\t" + str(Lc))


# CHECK NORMALIZATION

P1 = 0
P2 = 0
for i in range(N):
  P1 += TDMAT[i,i] / SUM1
  for j in range(N):
    P2 += TDMAT[i,j] / SUM2
print ("Check norm for probabilities:", P1,P2, "-- Both should be ~ 1.0")


