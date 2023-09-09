import numpy as np
from matplotlib import pyplot as plt

N = 42 # Dimension of Matrix

NStates = 20

for state in range( 1,NStates+1 ):
    lines = open(f'tmat_G_S{state}.txt','r').readlines()

    TDMat = np.zeros(( N,N ))
    for line in lines:
        t = line.split()
        ind1 = int( t[0] ) - 1
        ind2 = int( t[1] ) - 1
        TDMat[ind1,ind2] = float( t[2] )

    plt.imshow( TDMat, origin="lower")
    plt.savefig( f"tmat_G_S{state}.txt.jpg" )
