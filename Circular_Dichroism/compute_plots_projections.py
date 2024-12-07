import numpy as np
import subprocess as sp
import matplotlib
from matplotlib import pyplot as plt
from glob import glob

"""
    Intructions:
    1. Gaussian TDDFT with Additional Keywords: "IOp(6/8=3) IOp(9/40=5)"
    2. make_transition_density.slurm
    3. make_electric_dipole_density.slurm
    4. make_magnetic_dipole_density.slurm
    5. compute_plots_projections.py (Only required for reduced description of the dipole densities)
"""

def get_Globals():
    global NSTATES, PLOTS_DIR, dim_dict, CUBE_DIR
    CUBE_DIR="cube_files"

    PLOTS_DIR = "PLOTS"
    sp.call(f"mkdir -p {PLOTS_DIR}", shell=True)
    dim_dict = { 0:"X", 1:"Y", 2:"Z" }
    NSTATES = len(glob(f"{CUBE_DIR}/TRANSITION_DENSITY_0_*.cube"))
    print(NSTATES)

def get_DENSITY_DATA( state_j, state_k, DENS_TYPE ):
    """
    state_j (int)  : Electronic State Index
    state_k (int)  : Electronic State Index
    DENS_TYPE (str): "ELECTRIC" or "MAGNETIC" or "TRANSITION" 
                    -- Transition Dipole Density
                    -- Transition Magnetic Dipole Density
                    -- Transition Density
    """

    if ( DENS_TYPE == "ELECTRIC" ):
        file_name = "EL_DIP_DENSITY"
    elif( DENS_TYPE == "MAGNETIC" ):
        file_name = "MAG_DIP_DENSITY"
    elif( DENS_TYPE == "TRANSITION" ):
        file_name = "TRANSITION_DENSITY"
    else:
        print(f"Density type not valid: {DENS_TYPE}.")
        print("Must be 'TRANSITION' or 'ELECTRIC' or 'MAGNETIC'.")
        exit()

    # Get size from first TD cube file from QCHEM (not Gaussian)
    global NAtoms, NGrid, Nxyz, dLxyz, Lxyz, coords # I assume all cubes are the same size
    if ( DENS_TYPE == "TRANSITION" ):
        header = np.array([ np.array(j.split(),dtype=float) for j in open(f"{CUBE_DIR}/{file_name}_{state_j}_{state_k}.cube","r").readlines()[2:6] ])
    else:
        header = np.array([ np.array(j.split(),dtype=float) for j in open(f"{CUBE_DIR}/{file_name}_{state_j}_{state_k}_X.cube","r").readlines()[2:6] ])
    NAtoms = int(header[0,0])
    Nxyz   = (header[1:,0]).astype(int)
    NGrid  = Nxyz[0] * Nxyz[1] * Nxyz[2]
    dLxyz  = np.array([header[1,1],header[2,2],header[3,3] ]).astype(float)
    Lxyz   = np.array([ header[0,1], header[0,2], header[0,3] ]).astype(float)
    if ( Lxyz[0] < 0 ): 
        Lxyz  *= -1.000 # Switch Sign, Already Angstroms
    if ( Lxyz[0] > 0 ): 
        Lxyz  *= 0.529 # Convert from Bohr to Angstroms
        dLxyz *= 0.529 # Convert from Bohr to Angstroms
    Vol    = Lxyz[0] * Lxyz[1] * Lxyz[2]
    
    print (f'\tNAtoms      = {NAtoms}')
    print (f'\tNx Ny Nz    = {Nxyz[0]} {Nxyz[1]} {Nxyz[2]}')
    print (f'\tLx Ly Lz    = {Lxyz[0]} {Lxyz[1]} {Lxyz[2]} A')
    print (f'\tdLx dLy dLz = {dLxyz[0]} {dLxyz[1]} {dLxyz[2]} A')
    print (f'\tVolume      = {Vol} A^3')

    NStart = NAtoms + 6

    if ( DENS_TYPE == "TRANSITION" ):
        coords = np.array([ j for j in open(f"{CUBE_DIR}/{file_name}_{state_j}_{state_k}.cube","r").readlines()[6:NStart] ])
    else:
        coords = np.array([ j for j in open(f"{CUBE_DIR}/{file_name}_{state_j}_{state_k}_X.cube","r").readlines()[6:NStart] ])

    if ( DENS_TYPE == "TRANSITION" ):
        DENSITY = np.zeros(( Nxyz[0], Nxyz[1], Nxyz[2]  ))
    else:
        DENSITY = np.zeros(( Nxyz[0], Nxyz[1], Nxyz[2], 3  ))
    print(f"\tMemory size of transition density array in (MB, GB): ({round(DENSITY.size * DENSITY.itemsize * 10 ** -6,2)},{round(DENSITY.size * DENSITY.itemsize * 10 ** -9,2)})" )
    
    for dind,dim in enumerate(["X","Y","Z"]):
        print (f'\tReading Transition {DENS_TYPE} Dipole Density: {state_j}-{state_k} ({dim}).')
        temp = []
        try:
            if ( DENS_TYPE == "TRANSITION" ):
                lines = open(f"{CUBE_DIR}/{file_name}_{state_j}_{state_k}.cube",'r').readlines()[NStart:]
            else:
                lines = open(f"{CUBE_DIR}/{file_name}_{state_j}_{state_k}_{dim}.cube",'r').readlines()[NStart:]
        except FileNotFoundError:
            print (f'\t****** File "{CUBE_DIR}/{file_name}_{state_j}_{state_k}_{dim}.cube" not found. ******')
            continue
        for count, line in enumerate(lines):
            t = line.split('\n')[0].split()
            for j in range(len(t)):
                temp.append( float(t[j]) )
        if ( DENS_TYPE == "TRANSITION" ):
            DENSITY[:,:,:] = np.array( temp ).reshape(( Nxyz[0],Nxyz[1],Nxyz[2] ))
        else:
            DENSITY[:,:,:,dind] = np.array( temp ).reshape(( Nxyz[0],Nxyz[1],Nxyz[2] ))

        if ( np.allclose(DENSITY, np.zeros( DENSITY.shape )) ):
            print(state_j, state_k, "TD = ZEROS")
            exit()      
        if ( DENS_TYPE == "TRANSITION" ):
            break            

    return DENSITY

def write_3D_DATA( state_j, state_k, CD_0k, FILE_NAME ):
    f = open(f'{PLOTS_DIR}/{FILE_NAME}_{state_j}_{state_k}.cube','w')
    f.write(f"Circular Dichroism Density \n")
    f.write(f"Totally {NGrid} grid points\n")
    f.write(f"{NAtoms} {-Lxyz[0]/0.529} {-Lxyz[1]/0.529} {-Lxyz[2]/0.529}\n")
    f.write(f"{Nxyz[0]} {dLxyz[0]/0.529}  0.000000   0.000000\n")
    f.write(f"{Nxyz[1]} 0.000000   {dLxyz[1]/0.529} 0.000000\n")
    f.write(f"{Nxyz[2]} 0.000000   0.000000   {dLxyz[2]/0.529} \n")
    for at in range(len(coords)):
        f.write( coords[at] )
    for x in range(Nxyz[0]):
        #print(f'X = {x}')
        for y in range(Nxyz[1]):
            outArray = []
            for z in range(Nxyz[2]):
                outArray.append( CD_0k[x,y,z] )
                if ( len(outArray) % 6 == 0 or z == Nxyz[2]-1 ):
                    #outArray.append('\n')
                    f.write( " ".join(map( str, np.round(outArray,8) )) + "\n" )
                    outArray = []
    f.close()

def write_2D_DATA( state_j, state_k, DATA, FILENAME, d1, d2 ):
    np.savetxt(f"{PLOTS_DIR}/{FILENAME}_{state_j}_{state_k}.dat", DATA, fmt='%1.4f')

    RMIN1 = -Lxyz[d1]/2
    RMAX1 =  Lxyz[d1]/2
    RMIN2 = -Lxyz[d2]/2
    RMAX2 =  Lxyz[d2]/2
    norm = matplotlib.colors.TwoSlopeNorm(vcenter=0.0)
    plt.imshow( DATA, norm=norm, origin='lower', interpolation="spline16", cmap="bwr", extent=[RMIN1, RMAX1, RMIN2, RMAX2], aspect='equal' )
    plt.colorbar(pad=0.01)
    plt.xlabel(f"Position {dim_dict[d2]} ($\\AA$)", fontsize=15)
    plt.ylabel(f"Position {dim_dict[d1]} ($\\AA$)", fontsize=15)
    plt.tight_layout()
    plt.savefig(f"{PLOTS_DIR}/{FILENAME}_{state_j}_{state_k}.jpg", dpi=300)
    plt.clf()

def write_1D_DATA( state_j, state_k, DATA, FILENAME, d ):
    np.savetxt(f"{PLOTS_DIR}/{FILENAME}_{state_j}_{state_k}.dat", DATA, fmt='%1.4f')

    RMIN = -Lxyz[d]/2
    RMAX =  Lxyz[d]/2
    plt.plot( np.linspace(RMIN,RMAX,Nxyz[d]), DATA, color='black', linewidth=3 )
    plt.xlabel(f"Position {dim_dict[d]} ($\\AA$)", fontsize=15)
    plt.ylabel("Density", fontsize=15)
    plt.tight_layout()
    plt.savefig(f"{PLOTS_DIR}/{FILENAME}_{state_j}_{state_k}.jpg", dpi=300)
    plt.clf()

def run():
    for state_k in range( 1, NSTATES+1 ):
        TD_0k   = get_DENSITY_DATA( 0, state_k, "TRANSITION") # (Nx,Ny,Nz)
        EDIP_0k = get_DENSITY_DATA( 0, state_k, "ELECTRIC")   # (Nx,Ny,Nz,3)
        MDIP_0k = get_DENSITY_DATA( 0, state_k, "MAGNETIC")   # (Nx,Ny,Nz,3)
        
        print(f"\nTransition Properties (0->{state_k}):" )
        EL_MU_VAL  = np.einsum("xyzD->D",EDIP_0k) * np.prod(dLxyz/0.529)
        MAG_MU_VAL = np.einsum("xyzD->D",MDIP_0k) * np.prod(dLxyz/0.529)
        print("\tElectric Dipole (from EL_DIP):              \\vec{r} = %1.4f  %1.4f  %1.4f a.u." % (EL_MU_VAL[0], EL_MU_VAL[1], EL_MU_VAL[2]) )
        print("\tMagnetic Dipole (from MAG_DIP):    \\vec{r} x \\nabla = %1.4f  %1.4f  %1.4f a.u." % (MAG_MU_VAL[0], MAG_MU_VAL[1], MAG_MU_VAL[2] ) )

        ##### Read dipole data from gaussian output file (if in the present directory) #####
        try:
            G16_OUT = glob("*.out")[0]
            sp.call("grep -A %d 'electric dipole' %s | tail -n %d | awk '{print $2, $3, $4}' > temp.dat" % (NSTATES+1,G16_OUT,NSTATES), shell=True)
            G16_EDIP = np.loadtxt("temp.dat")
            sp.call("grep -A %d 'magnetic dipole' %s | tail -n %d | awk '{print $2, $3, $4}' > temp.dat" % (NSTATES+1,G16_OUT,NSTATES), shell=True)
            G16_MDIP = np.loadtxt("temp.dat")
            sp.call("rm temp.dat", shell=True)
            print("\tElectric Dipole (from G16 output):          \\vec{r} = %1.4f  %1.4f  %1.4f a.u." % (G16_EDIP[state_k-1,0],G16_EDIP[state_k-1,1],G16_EDIP[state_k-1,2]))
            print("\tMagnetic Dipole (from G16 output): \\vec{r} x \\nabla = %1.4f  %1.4f  %1.4f a.u." % (G16_MDIP[state_k-1,0],G16_MDIP[state_k-1,1],G16_MDIP[state_k-1,2]))
        except:
            pass

        ##### Compute the 1D Projections of the Various Densities (useful for 1D systems) #####
        # ELECTRIC dipole density (x_i)
        write_1D_DATA( 0, state_k, np.sum(TD_0k[:,:,:],axis=(1,2)), "TD_1D_X", 0)
        write_1D_DATA( 0, state_k, np.sum(TD_0k[:,:,:],axis=(0,2)), "TD_1D_Y", 1)
        write_1D_DATA( 0, state_k, np.sum(TD_0k[:,:,:],axis=(0,1)), "TD_1D_Z", 2)

        # ELECTRIC dipole density (x_i)
        write_1D_DATA( 0, state_k, np.sum(EDIP_0k[:,:,:,0],axis=(1,2)), "EL_DIP_1D_X", 0)
        write_1D_DATA( 0, state_k, np.sum(EDIP_0k[:,:,:,1],axis=(0,2)), "EL_DIP_1D_Y", 1)
        write_1D_DATA( 0, state_k, np.sum(EDIP_0k[:,:,:,2],axis=(0,1)), "EL_DIP_1D_Z", 2)

        # MAGNETIC dipole density (x_i)
        write_1D_DATA( 0, state_k, np.sum(MDIP_0k[:,:,:,0],axis=(1,2)), "MAG_DIP_1D_X", 0)
        write_1D_DATA( 0, state_k, np.sum(MDIP_0k[:,:,:,1],axis=(0,2)), "MAG_DIP_1D_Y", 1)
        write_1D_DATA( 0, state_k, np.sum(MDIP_0k[:,:,:,2],axis=(0,1)), "MAG_DIP_1D_Z", 2)


        ##### Compute various 2D Projections of the Various Densities (probably not useful) #####
        # TRANSITION density (x_i,x_j)
        write_2D_DATA( 0, state_k, np.einsum( "xyz->xy", TD_0k ), "TD_2D_xy", 0, 1 )
        write_2D_DATA( 0, state_k, np.einsum( "xyz->xz", TD_0k ), "TD_2D_xz", 0, 2 )
        write_2D_DATA( 0, state_k, np.einsum( "xyz->yz", TD_0k ), "TD_2D_yz", 1, 2 )

        # ELECTRIC dipole density (x_i,x_j)
        write_2D_DATA( 0, state_k, np.einsum( "xyz->xy", EDIP_0k[:,:,:,0] ), "EL_DIP_2D_X_xy", 0, 1 )
        write_2D_DATA( 0, state_k, np.einsum( "xyz->xz", EDIP_0k[:,:,:,0] ), "EL_DIP_2D_X_xz", 0, 2 )
        write_2D_DATA( 0, state_k, np.einsum( "xyz->yz", EDIP_0k[:,:,:,0] ), "EL_DIP_2D_X_yz", 1, 2 )
        write_2D_DATA( 0, state_k, np.einsum( "xyz->xy", EDIP_0k[:,:,:,1] ), "EL_DIP_2D_Y_xy", 0, 1 )
        write_2D_DATA( 0, state_k, np.einsum( "xyz->xz", EDIP_0k[:,:,:,1] ), "EL_DIP_2D_Y_xz", 0, 2 )
        write_2D_DATA( 0, state_k, np.einsum( "xyz->yz", EDIP_0k[:,:,:,1] ), "EL_DIP_2D_Y_yz", 1, 2 )
        write_2D_DATA( 0, state_k, np.einsum( "xyz->xy", EDIP_0k[:,:,:,2] ), "EL_DIP_2D_Z_xy", 0, 1 )
        write_2D_DATA( 0, state_k, np.einsum( "xyz->xz", EDIP_0k[:,:,:,2] ), "EL_DIP_2D_Z_xz", 0, 2 )
        write_2D_DATA( 0, state_k, np.einsum( "xyz->yz", EDIP_0k[:,:,:,2] ), "EL_DIP_2D_Z_yz", 1, 2 )


        # MAGNETIC dipole density (x_i,x_j)
        write_2D_DATA( 0, state_k, np.einsum( "xyz->xy", MDIP_0k[:,:,:,0] ), "MAG_DIP_2D_X_xy", 0, 1 )
        write_2D_DATA( 0, state_k, np.einsum( "xyz->xz", MDIP_0k[:,:,:,0] ), "MAG_DIP_2D_X_xz", 0, 2 )
        write_2D_DATA( 0, state_k, np.einsum( "xyz->yz", MDIP_0k[:,:,:,0] ), "MAG_DIP_2D_X_yz", 1, 2 )
        write_2D_DATA( 0, state_k, np.einsum( "xyz->xy", MDIP_0k[:,:,:,1] ), "MAG_DIP_2D_Y_xy", 0, 1 )
        write_2D_DATA( 0, state_k, np.einsum( "xyz->xz", MDIP_0k[:,:,:,1] ), "MAG_DIP_2D_Y_xz", 0, 2 )
        write_2D_DATA( 0, state_k, np.einsum( "xyz->yz", MDIP_0k[:,:,:,1] ), "MAG_DIP_2D_Y_yz", 1, 2 )
        write_2D_DATA( 0, state_k, np.einsum( "xyz->xy", MDIP_0k[:,:,:,2] ), "MAG_DIP_2D_Z_xy", 0, 1 )
        write_2D_DATA( 0, state_k, np.einsum( "xyz->xz", MDIP_0k[:,:,:,2] ), "MAG_DIP_2D_Z_xz", 0, 2 )
        write_2D_DATA( 0, state_k, np.einsum( "xyz->yz", MDIP_0k[:,:,:,2] ), "MAG_DIP_2D_Z_yz", 1, 2 )

        
        
        # Probably can stop here...Below is just for fun and debugging
        continue



        # What is Rotary Strength Density (or CD) Density
        """
        \\vec{MU}_EL  = \\int~dx dy dz  \\vec{x, y, z} \\rho_{TD}(x, y, z)
        \\vec{MU}_MAG = \\int~dx'dy'dz' \\vec{x',y',z'} \\times \\nabla \\rho_{TD}(x',y',z')
        CD_0k         = \\vec{MU}_EL \\cdot \\vec{MU}_MAG

        \\vec{MU}_EL(x, y, z)   = \\vec{x, y, z} \\rho_{TD}(x, y, z)
        \\vec{MU}_MAG(x',y',z') = \\vec{x',y',z'} x \\nabla \\rho_{TD}(x',y',z')
        CD_0k_2D(r,r') = \\vec{MU}_EL(x, y, z) * \\vec{MU}_MAG(x',y',z')
        CD_0k          = \\int~dx dy dz \\int~dx'dy'dz' \\vec{MU}_EL(x, y, z) * \\vec{MU}_MAG(x',y',z')
        """
        ### Total Rotary Strength Density is a 6D Object: R(r,r')
        ### Best we can do here is 2D projections. Not sure if this is even useful.
        ### R(x,x') = \int dy dz dy' dz'\vec{\mu}_{EL}(x,y,z) \cdot \vec{\mu}_{MAG}(x',y',z')
        ### There are 27 possible 2D projections of the 6D object R(r,r'). 
        ### I only did the diagonal ones: Xxx, Yyy, Zzz
        EDIP_0k_X = np.einsum( "xyzD->xD", EDIP_0k )
        EDIP_0k_Y = np.einsum( "xyzD->yD", EDIP_0k )
        EDIP_0k_Z = np.einsum( "xyzD->zD", EDIP_0k )
        MDIP_0k_X = np.einsum( "xyzD->xD", MDIP_0k )
        MDIP_0k_Y = np.einsum( "xyzD->yD", MDIP_0k )
        MDIP_0k_Z = np.einsum( "xyzD->zD", MDIP_0k )
        write_2D_DATA( 0, state_k, np.einsum( "aD,xD->ax", EDIP_0k_X, MDIP_0k_X ), "CD_2D_X_xx", 0, 0 )
        write_2D_DATA( 0, state_k, np.einsum( "bD,yD->by", EDIP_0k_Y, MDIP_0k_Y ), "CD_2D_Y_yy", 1, 1 )
        write_2D_DATA( 0, state_k, np.einsum( "cD,zD->cz", EDIP_0k_Z, MDIP_0k_Z ), "CD_2D_Z_zz", 2, 2 )
        



        ##### FOR DEBUGGING -- REWRITE THE DATA WE ALREADY HAVE #####
        ## Write TRANSITION density to cube file ##
        # write_3D_DATA( 0, state_k, TD_0k[:,:,:], "TD_3D" )

        ## Write ELECTRIC dipole density to cube file ##
        # write_3D_DATA( 0, state_k, EDIP_0k[:,:,:,0], "EL_DIP_3D_X" )
        # write_3D_DATA( 0, state_k, EDIP_0k[:,:,:,1], "EL_DIP_3D_Y" )
        # write_3D_DATA( 0, state_k, EDIP_0k[:,:,:,2], "EL_DIP_3D_Z" )

        ## Write MAGNETIC dipole density to cube file ##
        # write_3D_DATA( 0, state_k, MDIP_0k[:,:,:,0], "MAG_DIP_3D_X" )
        # write_3D_DATA( 0, state_k, MDIP_0k[:,:,:,1], "MAG_DIP_3D_Y" )
        # write_3D_DATA( 0, state_k, MDIP_0k[:,:,:,2], "MAG_DIP_3D_Z" )




def main():
    get_Globals()
    run()


if ( __name__ == "__main__" ):
    main()