import numpy as np
from matplotlib import pyplot as plt



def get_Globals():
    global NROOTS, TDMat_DIR, SWCNT_AXIS_DIR
    NROOTS = 20
    TDMat_DIR = "./" # Location of the TransDens files





def get_TD_Data():

    #print ("\tStarting to Read TD Files.")

    # Get size from first TD file
    global NAtoms, NGrid, Nxyz, dLxyz, Lxyz, coords
    header = np.array([ j.split() for j in open(f"{TDMat_DIR}/trans-0_1.cube","r").readlines()[1:6] ])
    NAtoms = int(header[1,0])
    NGrid  = int( header[0,1] )
    Nxyz   = (header[2:,0]).astype(int)
    dLxyz  = np.array([header[2,1],header[3,2],header[4,3] ]).astype(float)
    #Lxyz   = np.array([ Nxyz[0]*dLxyz[0], Nxyz[1]*dLxyz[1], Nxyz[2]*dLxyz[2] ])
    Lxyz   = np.array([ header[1,1], header[1,2], header[1,3] ]).astype(float)
    if ( Lxyz[0] < 0 ): Lxyz *= -1.000 # Already Angstroms
    elif ( Lxyz[0] > 0 ): Lxyz *= 0.529 # Convert from Bohr to Angstroms
    Vol    = Lxyz[0] * Lxyz[1] * Lxyz[2]
    

    print (f'\tNAtoms      = {NAtoms}')
    print (f'\tTD Grid     = {NGrid}')
    print (f'\tNx Ny Nz    = {Nxyz[0]} {Nxyz[1]} {Nxyz[2]}')
    print (f'\tLx Ly Lz    = {Lxyz[0]} {Lxyz[1]} {Lxyz[2]} A')
    print (f'\tdLx dLy dLz = {dLxyz[0]} {dLxyz[1]} {dLxyz[2]} A')
    print (f'\tVolume      = {round(Vol,2)} A^3')

    NStart = NAtoms + 6

    coords = np.array([ j.strip("\n").split()[1:] for j in open(f"{TDMat_DIR}/trans-0_1.cube","r").readlines()[6:NStart] ]).astype(float)
    coords *= 0.529 # Coords are in Bohr still, even when voxels are in Angstrom... Is this always true ?

    TD = np.zeros(( NROOTS, Nxyz[0], Nxyz[1], Nxyz[2]  ))
    print(f"\tMemory size of transition density array in (MB, GB): ({round(TD.size * TD.itemsize * 10 ** -6,2)},{round(TD.size * TD.itemsize * 10 ** -9,2)})" )
    for state_j in range(NROOTS):
        #print (f'\tReading Transition Density {mat+1} of {NROOTS}.')
        temp = []
        try:
            lines = open(f"{TDMat_DIR}/trans-0_{state_j+1}.cube",'r').readlines()[NStart:]
        except FileNotFoundError:
            print (f'\t****** File "trans-0_{state_j+1}.cube" not found. Skipping this matrix element. ******')
            continue
        for count, line in enumerate(lines):
            t = line.split('\n')[0].split()
            for j in range(len(t)):
                temp.append( float(t[j]) )
        TD[state_j,:,:,:] = np.array( temp ).reshape(( Nxyz[0],Nxyz[1],Nxyz[2] ))
        #compute_Electrostatic_Moments( TD, state_j )                   
                    
    return TD

def compute_Ld(TD):
    Ld = np.zeros(( NROOTS, 2 )) # Angstrom, % SWCNT axis
    global SWCNT_AXIS_DIR
    SWCNT_AXIS_DIR = np.argmax( Nxyz ) # THIS IS HIGHLY SYSTEM DEPENDENT.
    SWCNT_L = np.max( coords[:,SWCNT_AXIS_DIR+1] ) - np.min( coords[:,SWCNT_AXIS_DIR+1] )
    print( "\t\tSWCNT Length (A):", SWCNT_L )
    print(f"\t\tCHOOSING MAIN AXIS = {SWCNT_AXIS_DIR}")
    TD_NORMED = np.zeros(( NROOTS, Nxyz[SWCNT_AXIS_DIR] ))
    for state in range( NROOTS ):

        # Get probability
        axis               = (1,2)*(SWCNT_AXIS_DIR==0)+(2,0)*((SWCNT_AXIS_DIR==1))+(0,1)*((SWCNT_AXIS_DIR==2))
        TD_NORMED[state,:] = np.sum( np.abs(TD[state,:,:,:]), axis=axis )
        NORM               = np.sum( np.abs(TD_NORMED[state,:]) ) # Total Population
        TD_NORMED[state,:] = np.abs(TD_NORMED[state,:]) / NORM

        # Test Compute Ld
        #Ld[state,:] = 1 / np.sum( np.array([1,0,0,0,0,0,0,0,0,0]) ) # TEST: Ld --> 1
        #Ld[state,:] = 1 / np.sum( (np.ones(11)/11)**2 ) # TEST: Ld --> 10
        # Compute Ld
        Ld[state,0] = 1 / np.sum( TD_NORMED[state,:]**2 )
        Ld[state,0] *= dLxyz[SWCNT_AXIS_DIR] * 0.529 # [0,1] Boxes --> Bohr --> Angstroms
        Ld[state,1]  = Ld[state,0] * dLxyz[SWCNT_AXIS_DIR] / SWCNT_L * 100 # [0,1] --> Angstroms --> % SWCNT Axis
        
        print( f"State {state+1}; Ld (A/%) =", round(Ld[state,0],2), "/", round(Ld[state,1],2) )

    # Save Ld
    np.savetxt("Ld.dat", Ld, fmt="%2.5f", header="Ld (ANG)    Ld(% SWCNT Axis)")

    # Save normalized TD
    np.savetxt("TD_NORM.dat", TD_NORMED[:,:].T, fmt="%2.5f")

    return TD_NORMED, Ld

def make_plots(TD_NORMED, Ld):
    LGrid = np.arange(Nxyz[SWCNT_AXIS_DIR])*dLxyz[SWCNT_AXIS_DIR]*0.529
    MAX_TD = np.max( TD_NORMED )
    for state in range( 10 ):
        plt.plot( LGrid, state*MAX_TD*np.ones(len(LGrid)), c='black', lw=1 )
        plt.plot( LGrid, state*MAX_TD+TD_NORMED[state,:], lw=3, label=f"S{state+1}" )
        plt.text( LGrid[-1]+1,(state+0.2)*MAX_TD,f"L$_d$ = {int(np.round(Ld[state,0],0))} $\AA$" )

    #plt.legend()
    plt.tick_params(left = False, labelleft = False, top = True)
    plt.xlim(LGrid[0],LGrid[-1])
    plt.ylim(0)
    plt.xlabel("SWCNT Axis Position, X ($\AA$)", fontsize=15)
    plt.ylabel("Transition Density, |$\\xi$ (X)| ($\AA$)", fontsize=15)
    plt.tight_layout()
    plt.savefig("TransDens_1D.jpg",dpi=600)

def get_moments( TD_NORMED ):
    """
    a_th Moment: M(a) = int[ X^a * TD(X) ] / int[ TD(X) ] = int[ X^a * TD_NORMED(X) ]
    """
    #print("Starting calculation of moments:")
    LGrid = np.arange(Nxyz[SWCNT_AXIS_DIR])*dLxyz[SWCNT_AXIS_DIR]
    M = np.zeros(( 3, NROOTS ))
    
    M[0,:] = np.einsum("b,sb->s", LGrid[:]*0.0 + 1 , TD_NORMED[:,:] )
    M[1,:] = np.einsum("b,sb->s", LGrid[:]    , TD_NORMED[:,:] )
    M[2,:] = np.einsum("b,sb->s", LGrid[:]**2 , TD_NORMED[:,:] )

    #print( "M(0):", np.round(M[0,:],0) )
    #print( "M(1):", np.round(M[1,:],0) )
    #print( "M(2) - M(1)^2:", np.round( M[2,:] - M[1,:]**2,4) )

    np.savetxt( "Moments.dat", M[:,:], fmt="%2.2f" )

if ( __name__ == "__main__" ):
    get_Globals()
    TD = get_TD_Data()
    TD_NORMED, Ld = compute_Ld(TD)
    get_moments(TD_NORMED)
    make_plots(TD_NORMED,Ld)
