# The purpose of this code is to read in tdmat_G_S1.txt
#  or an equivalent file produced by Multiwfn and extract 
#  two "slices" of the matrix for plotting: Ld and Lc

# SYNTAX: python3 plot_slice_TDMAT.py tdmat_G_S1.txt (3-Column File)

import sys
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

file01 = open("slice_width.txt","r")
dL = float(file01.readlines()[0])
file01.close()

name = sys.argv[1]
file01 = open(name,'r')
N = 0
for count, line in enumerate(file01):
  N = count
N = int(N ** (1/2))+1 # Dimension of Matrix
print ("Dimension of the transition density matrix:", N)
file01.close()

file01 = open(name,'r')
TDMAT = np.zeros((N,N))
for count, line in enumerate(file01):
  s = line.split()
  n = int(s[0])-1
  m = int(s[1])-1
  val = float(s[2])
  TDMAT[n,m] = val

# Make output file for diagonal (Ld), anti-diagonal (Lc),
#  and vertical (Lc?), all through the midpoint.

Ld = np.zeros((N))
Lc = Ld
vert = Ld
horiz = Ld

Ld = np.array( [ TDMAT[j,j] for j in range(N) ] )
Lc = np.array( [ TDMAT[N - j - 1,j] for j in range(N) ] ) # Check this later


xvals = np.arange(N) * dL / 10 # Old X-Points
xvals_new = np.linspace(0,np.max(xvals),500) # New x-points

# Do Ld
Ld_interp = interp1d(xvals, Ld, kind='cubic') # Interpolate
Ld_smooth = savgol_filter(Ld_interp(xvals_new), 51, 3) # Smooth
# Do Lc
Lc_interp = interp1d(xvals, Lc, kind='cubic') # Interpolate
Lc_smooth = savgol_filter(Lc_interp(xvals_new), 51, 3) # Smooth

# Plot Ld
plt.plot( xvals, Ld, 'o', label="Raw")
plt.plot( xvals_new, Ld_interp(xvals_new), "--", label="Interp")
plt.plot( xvals_new, Ld_smooth, "-", label="Smooth" )
plt.legend()
plt.savefig(name+"_Ld.jpg")
plt.clf()

# Plot Lc
plt.plot( xvals, Ld, 'o', label="Raw")
plt.plot( xvals_new, Ld_interp(xvals_new), "--", label="Interp")
plt.plot( xvals_new, Ld_smooth, "-", label="Smooth" )
plt.legend()
plt.savefig(name+"_Lc.jpg")
plt.clf()

# Plot Ld, Lc
plt.plot( xvals_new, Ld_smooth, "-", label="Ld")
plt.plot( xvals_new, Lc_smooth, "--", label="Lc" )
plt.legend()
plt.savefig(name+"_Ld_Lc.jpg")
plt.clf()

np.savetxt( name + "_diagonal_Ld.txt_raw", np.array( [xvals, Ld] ).T )
np.savetxt( name + "_diagonal_Lc.txt_raw", np.array( [xvals, Lc] ).T )
np.savetxt( name + "_diagonal_Ld.txt_smoothed", np.array( [xvals_new, Ld_smooth] ).T )
np.savetxt( name + "_diagonal_Lc.txt_smoothed", np.array( [xvals_new, Lc_smooth] ).T )

# Plot the full matrix
plt.contourf( TDMAT )
plt.xlabel("Electron (A)")
plt.ylabel("Hole (A)")
plt.title(name)
plt.colorbar()
plt.savefig( name + "_TDMAT.jpg")


"""

for i in range(N):
  for j in range(N):
    if (i == j): # Diagonal Line
      #Ld[i] = TDMAT[i,j]
      X = (j+1) * dL / 10 # Convert index to nm
      file01.write(str(X) + '\t' + str(TDMAT[i,j]) + '\n')
    if (i == N - j - 1): # Anti-diagonal Line
      #Lc[i] = TDMAT[i,j]
      X = (i+1)*dL/10
      file02.write(str(X) + '\t' + str(TDMAT[i,j]) + '\n')
    if (i == int(N/2)): # Vertical Line
      #vert[i] = TDMAT[i,j]
      X = (j+1)*dL/10
      file03.write(str(X) + '\t' + str(TDMAT[i,j]) + '\n')
    if (j == int(N/2)): # Horizontal  Line
      #horiz[i] = TDMAT[i,j]
      X = (i+1)*dL/10
      file04.write(str(X) + '\t' + str(TDMAT[i,j]) + '\n')


file01.close()
file01.close()
file03.close()
file04.close()
"""

