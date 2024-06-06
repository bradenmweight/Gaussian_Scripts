import numpy as np
from matplotlib import pyplot as plt

def get_Globals():
    global NGEOMS, TD_DIRS
    NGEOMS = 600
    TD_DIRS = [ f"TDDFT/STEP_{j}/PLOTS_DATA/" for j in range( NGEOMS ) ]

    global EMIN, EMAX, NPTS, SIG, EGRID
    EMIN  = 4.0 # eV
    EMAX  = 8.0 # eV
    SIG   = 0.025 # eV
    NPTS  = 1000
    EGRID = np.linspace( EMIN, EMAX, NPTS )

def get_data():
    NSTATES = len( np.loadtxt( f"{TD_DIRS[0]}/ADIABATIC_ENERGIES_RPA.dat" ) )

    ENERGY  = np.zeros( (NGEOMS,NSTATES) )
    OSC     = np.zeros( (NGEOMS,NSTATES) )

    for step in range( NGEOMS ):
        ENERGY[step,:] = np.loadtxt( f"{TD_DIRS[step]}/ADIABATIC_ENERGIES_RPA.dat" )[:NSTATES]
        DIPOLE         = np.load( f"{TD_DIRS[step]}/DIPOLE_RPA.dat.npy" )[0,:NSTATES,:] # g -> e, xyz
        dE             = ENERGY[step,:] - ENERGY[step,0]
        OSC[step,:]   += (2/3) * dE * np.einsum("Jd,Jd->J", DIPOLE, DIPOLE )
    
    ENERGY *= 27.2114 # eV

    return ENERGY, OSC

def plot_spectrum( ENERGY, OSC ):

    # Shift to transition energy
    dE  = ENERGY[:,1:] - ENERGY[:,0][:,None]
    OSC = OSC[:,1:]

    ABS = np.zeros( NPTS )
    for pt,E in enumerate( EGRID ):
        ABS[pt] = np.sum( OSC[:,:] * np.exp(-(E-dE[:,:])**2/2/SIG**2) )
    
    plt.plot( EGRID, ABS/NGEOMS )
    for step in range( NGEOMS ):
        plt.scatter( dE[step,:], OSC[step,:], s=1, c="red" )
        # markerline, stemlines, baseline = plt.stem( dE[step,:], OSC[step,:], linefmt="blue", markerfmt="o" )
        # markerline.set_markerfacecolor('none')
        # markerline.set_markeredgecolor('blue')
        # markerline.set_markersize(1)
        # markerline.set_markeredgewidth(0.1)

    plt.xlabel("Energy (eV)", fontsize=15)
    plt.ylabel("Absorption (Osc. Str.)", fontsize=15)
    plt.xlim(EMIN,EMAX)
    plt.ylim(0)
    plt.tight_layout()
    plt.savefig( "%s/SPEC_ENSEMBLE.jpg" % (TD_DIRS[0].split("/")[0]), dpi=300  )
    plt.clf()

def plot_Sn( ENERGY, n=1 ):
    plt.plot( np.arange(NGEOMS), ENERGY[:,0] - ENERGY[0,0], label="$S_0$" )
    for state in range( 1, n+1 ):
        plt.plot( np.arange(NGEOMS), ENERGY[:,state] - ENERGY[0,0], label="$S_%s$" % state )
    plt.xlim(0,NGEOMS)
    #plt.ylim(0)
    plt.xlabel("Reaction Step", fontsize=15)
    plt.ylabel("Energy (eV)", fontsize=15)
    plt.legend()
    plt.tight_layout()
    plt.savefig( "%s/S0_S1.jpg" % (TD_DIRS[0].split("/")[0]), dpi=300  )
    plt.clf()

def main():
    get_Globals()
    ENERGY, OSC = get_data()
    plot_spectrum( ENERGY, OSC )
    plot_Sn( ENERGY, n=2 )

if ( __name__ == "__main__" ):
    main()

