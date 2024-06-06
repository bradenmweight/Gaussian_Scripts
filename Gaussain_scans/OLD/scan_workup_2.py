# The purpose of this script is to read in TD_*/geometry.out
#   files and extract spectra (waterfall plot), Ei vs SP,
#   and Osc Str. vs SP#

# ~ Braden M. Weight, July 11, 2020  ~ #

# PYTHON NAME: scan_workup_2.py

import numpy as np
import subprocess as sp
import os
from matplotlib import pyplot as plt

NSteps = 20 # TD_0, .., TD_{NSteps}

# Check if we are doing spectra (waterfall plot and 2D E and Osc.Str. vs Scan#)
doSpectra = True
QLogZ = False # Take LOG of OSC.STR. before plotting -- See minute features -- Usually turn off
spectraFile = "~/Gaussain_scripts/Excited_States/Gen_spectra_gaus16.py"

### CHANGE THESE PARAMETERS IF MAKING SPECTRA
sigma = 0.05 # Gaussian STD in meV
emin = 2
emax = 8
npts = 2000
###

# These are not ready yet...
doTDMandNTOs = True # Submit transition density and NTOs?
makeTD_SUBMITFILE = "~/Gaussain_scripts/Transition_Density/make_rho_ij_low_density.slurm" # SBATCH FILE
makeNTOS_SUBMITFILE = "~/Gaussain_scripts/NTOs/submit.gaussian" # SBATCH FILE


##### MAIN CODE #####

# Read in data from TD-DFT output
dirs = [ f"TD_{j}" for j in range(5,NSteps+1) ]

def getSpec(data):

    dE = (emax - emin) / npts
    spec = np.zeros(( len(data), npts, 2 ))
    for j in range( len(data) ):
        for i in range(npts):
            x = emin + i*dE
            y = 0
            spec[j,i,0] = x
            spec[j,i,1] = np.sum( data[j,:,1] * np.exp(-(x-data[j,:,0])**2/2/sigma**2) ) # Sum over excitations
    return spec


if ( doSpectra ):
    for d in range( len(dirs) ):
        print (f"Dir: {dirs[d]}")
        os.chdir(f"{dirs[d]}/")
        sp.call(f"python3 {spectraFile}", shell=True)
        if ( d == 0):
            global data
            tmp = np.loadtxt("Energy_Osc.dat")
            data = np.zeros(( NSteps+1, len(tmp), 2 )) # Scan Points + 1, NStates from TDDFT, energy/osc
            data[d] = tmp
        else:
            data[d] = np.loadtxt("Energy_Osc.dat")
        os.chdir("../")

    ### Make contourf plot for spectra
    spec = getSpec(data) # NGeoms, NPTS, Energy/Osc.Str.
    X,Y = np.meshgrid( spec[0,:,0], np.arange( len(spec) ) )
    if ( QLogZ ):
        plt.contourf( X,Y, np.log( np.abs( np.round( spec[:,:,1] ,4) ) ) )
        cbar = plt.colorbar()
        cbar.ax.set_title("Log( Osc. Str.)")
    else:
        plt.contourf( X,Y, np.round(spec[:,:,1],2) )
        cbar = plt.colorbar()
        cbar.ax.set_title("Oscillator Strength")
    plt.xlabel("Photon Energy (eV)", fontsize=15)
    plt.ylabel("Scan Step #", fontsize=15)
    plt.savefig("contourf.jpg")
    plt.clf()

    ### Make 2, 2D figures: Energy vs. Scan# and Osc. Str. vs. Scan#
    fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1, sharex=True)
    fig.subplots_adjust(hspace=0)
    for state in range( 5 ):
        ax1.plot( np.arange( len(data) ), data[:,state,0],"-o" )
        ax2.plot( np.arange( len(data) ), data[:,state,1],"-o", label=f"S{state+1}" )
    ax1.set_ylabel("Photon Energy (eV)", fontsize=12)
    ax2.set_ylabel("Oscillator Strength", fontsize=12)
    ax2.set_xlabel("Scan Step #", fontsize=15)
    ax2.legend()
    
    
    plt.savefig("E_Osc_vs_ScanPoint.jpg")
    plt.clf()

if ( doTDMandNTOs ):

    QDelete = False

    for SP in range( len(dirs) ):

        print (f"Dir: {dirs[SP]} ({SP} of {len(dirs)})")

        # Move into the TD-DFT directory
        os.chdir( dirs[SP] )

        # Make transition density matrix directory and move into it
        if ( os.path.exists("TDMat") and QDelete == False ):
            print (f"\tFound dirs: {dirs[SP]}/TDMat/. Shoul I remove all? [y,n]")
            QDeleteA = input()
            if ( QDeleteA == 'y' ):
                print ("\tRemoving directories: TDMat")
                sp.call("rm -r TDMat/",shell=True)
                QDelete = True
            else:
                print (f"\tQuitting. User said {QDeleteA}.")
                exit()
        elif ( os.path.exists("TDMat") and QDelete ):
            sp.call("rm -r TDMat/",shell=True)

        sp.call("mkdir TDMat/",shell=True)
        os.chdir( "TDMat/" )

        # Copy files, modify active state, and submit
        sp.call(f"cp {makeTD_SUBMITFILE} .",shell=True)
        if ( os.path.isfile("../geometry.fchk") ):
            sp.call(f"sbatch {makeTD_SUBMITFILE.split('/')[-1]}",shell=True)
        else:
            print  ("'geometry.fchk' does not exist. Submit TDMat files manually after making this.")

        os.chdir("../")

        # Make NTO directory and move into it
        if ( os.path.exists("NTOs") ):
            sp.call("rm -r NTOs/",shell=True)

        sp.call("mkdir NTOs/",shell=True)
        os.chdir( "NTOs/" )

        # Copy files, modify active state, and submit
        sp.call(f"cp {makeNTOS_SUBMITFILE} .",shell=True)
        sp.call(f"sbatch {makeNTOS_SUBMITFILE.split('/')[-1]}",shell=True)

        os.chdir("../../")

