# The purpose of this code is to slice the tube along its axis
#  and output a file to be read by Multiwfn to calculate the
#  fragment transition density of the system for some excited 
#  state.

import numpy as np

dL = 0.5 # slice width in Angstroms

file05 = open("slice_width.txt","w")
file05.write(str(dL))
file05.close()

# READ IN ATOMS
file01 = open("geometry.com_FULL","r")
atoms = []
flag = False
for count, line in enumerate(file01):
  s = line.split()
  if (len(s) == 2 and s[0] == "0" and s[1] == "1"):
    flag = True
  elif (len(s) == 2 and s[0] == "0" and s[1] == "2"):
    flag = True
  elif (len(s) == 2 and s[0] == "0" and s[1] == "3"):
    flag = True
  if (flag == True and len(s) == 4):
    x = float(s[1])
    y = float(s[2])
    z = float(s[3])
    atoms.append([x,y,z])
file01.close()
print ("Number of Atoms:", len(atoms))

# FIND TUBE DIRECTION
atoms = np.array(atoms)
xmax = np.amax(atoms, axis=0)[0]
ymax = np.amax(atoms, axis=0)[1]
zmax = np.amax(atoms, axis=0)[2]
xmin = np.amin(atoms, axis=0)[0]
ymin = np.amin(atoms, axis=0)[1]
zmin = np.amin(atoms, axis=0)[2]

Lx = xmax - xmin
Ly = ymax - ymin
Lz = zmax - zmin
L = [Lx,Ly,Lz]

tube_axis = []
if (Lx > Ly and Lx > Lz):
  tube_axis = [1,0,0]
if (Ly > Lx and Ly > Lz):
    tube_axis = [0,1,0]
if (Lz > Lx and Lz > Ly):
    tube_axis = [0,0,1]

print ("Tube axis is along", tube_axis)
print ("and of length", np.dot(L,tube_axis), "A")

# Find number of slices N for interval dL
N = int(np.dot(L,tube_axis) / dL)
print ("There will be N =", N, "slices")
print ("for interval dL =", dL, "A")

# Bin atoms by slice
T0 = np.dot(np.amin(atoms, axis=0),tube_axis) # Tube start
bins = [ [] for i in range(N) ]
for b in range(N):
  for n in range(len(atoms)):
    x = np.dot(atoms[n],tube_axis) # Atom coordinate
    if ( x > T0 + dL*b and x < T0 + dL*(b+1)):
      bins[b].append(n)
#for b in bins:
  #print (len(b), b)

# Export to file readable by Multiwfn
file01 = open("fragdef.txt","w")
for b in range(N):
  if (len(bins[b]) == 0):
    continue
  if (b != 0 and b != N):
    file01.write("\n")
  for n in range(len(bins[b])):
    file01.write(str(bins[b][n]))
    if ( n != len(bins[b])-1):
      file01.write(", ")




























