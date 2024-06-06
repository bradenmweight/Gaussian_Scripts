import numpy as np
import subprocess as sp
from matplotlib import pyplot as plt

"""
This script computes the rotary strength density from the 
    transition dipole and transition magnetic dipole densities
    as output from MultiWfn (V3.7).
    
    ~~~ Braden M. Weight -- May 9, 2024 ~~~

    Intructions:
    1. Gaussian TDDFT with Additional Keywords: "IOp(6/8=3) IOp(9/40=5)"
    2. make_dipole_density.slurm
    3. make_magnetic_dipole_density.slurm
    4. compute_CD_Density.py
"""

def get_Globals():
    global NSTATES
    NSTATES = 10



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
        file_name = "DIP_DENSITY"
    elif( DENS_TYPE == "MAGNETIC" ):
        file_name = "MAG_DIP_DENSITY"
    elif( DENS_TYPE == "TRANSITION" ):
        file_name = "TRANSITION_DENSITY"
    else:
        print(f"Density type not valid: {DENS_TYPE}.")
        print("Must be either 'ELECTRIC' or 'MAGNETIC'.")
        exit()

    # Get size from first TD cube file from QCHEM (not Gaussian)
    global NAtoms, NGrid, Nxyz, dLxyz, Lxyz, coords
    if ( DENS_TYPE == "TRANSITION" ):
        header = np.array([ np.array(j.split(),dtype=float) for j in open(f"{file_name}_{state_j}_{state_k}.cube","r").readlines()[2:6] ])
    else:
        header = np.array([ np.array(j.split(),dtype=float) for j in open(f"{file_name}_{state_j}_{state_k}_X.cube","r").readlines()[2:6] ])
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
        coords = np.array([ j for j in open(f"{file_name}_{state_j}_{state_k}.cube","r").readlines()[6:NStart] ])
    else:
        coords = np.array([ j for j in open(f"{file_name}_{state_j}_{state_k}_X.cube","r").readlines()[6:NStart] ])

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
                lines = open(f"{file_name}_{state_j}_{state_k}.cube",'r').readlines()[NStart:]
            else:
                lines = open(f"{file_name}_{state_j}_{state_k}_{dim}.cube",'r').readlines()[NStart:]
        except FileNotFoundError:
            print (f'\t****** File "{file_name}_{state_j}_{state_k}_{dim}.cube" not found. ******')
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
    f = open(f'{FILE_NAME}_{state_j}_{state_k}.cube','w')
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

def write_2D_DATA( state_j, state_k, DATA, FILENAME, d ):
    np.savetxt(f"{FILENAME}_{state_j}_{state_k}.dat", DATA, fmt='%1.4f')

    RMIN = -Lxyz[d]/2
    RMAX =  Lxyz[d]/2
    plt.imshow( DATA, origin='lower', interpolation="spline16", cmap="afmhot_r", extent=[RMIN, RMAX, RMIN, RMAX] )
    plt.colorbar(pad=0.01)
    plt.savefig(f"{FILENAME}_{state_j}_{state_k}.png", dpi=300)
    plt.clf()

def write_1D_DATA( state_j, state_k, DATA, FILENAME, d ):
    np.savetxt(f"{FILENAME}_{state_j}_{state_k}.dat", DATA, fmt='%1.4f')

    RMIN = -Lxyz[d]/2
    RMAX =  Lxyz[d]/2
    plt.plot( np.linspace(RMIN,RMAX,Nxyz[d]), DATA )
    plt.savefig(f"{FILENAME}_{state_j}_{state_k}.png", dpi=300)
    plt.clf()

def compute_CD_Density():
    #for state_k in range( 1, NSTATES+1 ):
    #for state_k in [1,3,8]:
    #for state_k in [1,3]:
    for state_k in [8]:
        TD_0k   = get_DENSITY_DATA( 0, state_k, "TRANSITION") # (Nx,Ny,Nz)
        EDIP_0k = get_DENSITY_DATA( 0, state_k, "ELECTRIC") # (Nx,Ny,Nz,3)
        MDIP_0k = get_DENSITY_DATA( 0, state_k, "MAGNETIC") # (Nx,Ny,Nz,3)
        
        xGRID = np.linspace( -Lxyz[0]/2, Lxyz[0]/2, Nxyz[0] ) * 0.529
        yGRID = np.linspace( -Lxyz[1]/2, Lxyz[1]/2, Nxyz[1] ) * 0.529
        zGRID = np.linspace( -Lxyz[2]/2, Lxyz[2]/2, Nxyz[2] ) * 0.529

        print(f"\nTransition Properties (0->{state_k}):" )
        MU_X_EL = np.einsum( "xyz,x->", TD_0k, -xGRID ) * dLxyz[0] * 0.529
        MU_Y_EL = np.einsum( "xyz,y->", TD_0k, -yGRID ) * dLxyz[1] * 0.529
        MU_Z_EL = np.einsum( "xyz,z->", TD_0k, -zGRID ) * dLxyz[2] * 0.529
        print("\tElectric Dipole (TD):                  \\vec{r} = %1.3f  %1.3f  %1.3f" % (MU_X_EL, MU_Y_EL, MU_Z_EL) )
        MU_X_EL = np.sum(EDIP_0k[:,:,:,0]) * 0.529**4 * dLxyz[0]
        MU_Y_EL = np.sum(EDIP_0k[:,:,:,1]) * 0.529**4 * dLxyz[1]
        MU_Z_EL = np.sum(EDIP_0k[:,:,:,2]) * 0.529**4 * dLxyz[2]
        print("\tElectric Dipole (EL DIP):              \\vec{r} = %1.3f  %1.3f  %1.3f" % (MU_X_EL, MU_Y_EL, MU_Z_EL) )
        MU_X_MAG = np.sum(MDIP_0k[:,:,:,0]) * 0.529**4 * dLxyz[0]
        MU_Y_MAG = np.sum(MDIP_0k[:,:,:,1]) * 0.529**4 * dLxyz[1]
        MU_Z_MAG = np.sum(MDIP_0k[:,:,:,2]) * 0.529**4 * dLxyz[2]
        print("\tMagnetic Dipole (MAG DIP):    \\vec{r} x \\nabla = %1.3f  %1.3f  %1.3f" % (MU_X_MAG, MU_Y_MAG, MU_Z_MAG) )
        # R_X = MU_X_EL * MU_X_MAG
        # R_Y = MU_Y_EL * MU_Y_MAG
        # R_Z = MU_Z_EL * MU_Z_MAG
        # R   = np.sqrt( R_X**2 + R_Y**2 + R_Z**2 ) * (4.6596/27.2114)
        # print("\tRotary Strength (Length): \mu_EL \cdot \mu_MAG = %1.3f" % (R) )
        # print( 3.5 / R )
        # # We need erg-esu-cm/Gauss
        # length_cgs = 5.29  * 10**(-9)   # cm
        # E_cgs      = 4.36  * 10**(-11)  # ergs
        # Q_cgs      = 4.803 * 10**(-10)  # esu
        # mag_cgs    = 1.400 * 10**6      # Mhz/Gauss
        # print(  R * E_cgs * Q_cgs * length_cgs * mag_cgs / 10**(-40)  )


        # Write transition density to cube file
        write_3D_DATA( 0, state_k, TD_0k[:,:,:], "TD_3D" )
        write_1D_DATA( 0, state_k, np.sum(TD_0k[:,:,:],axis=(0,1)), "TD_1D_X", 0)
        write_1D_DATA( 0, state_k, np.sum(TD_0k[:,:,:],axis=(0,2)), "TD_1D_Y", 1)
        write_1D_DATA( 0, state_k, np.sum(TD_0k[:,:,:],axis=(0,1)), "TD_1D_Z", 2)
        
        # Write ELECTRIC dipole density to cube file
        write_3D_DATA( 0, state_k, EDIP_0k[:,:,:,0], "EL_DIP_3D_X" )
        write_3D_DATA( 0, state_k, EDIP_0k[:,:,:,1], "EL_DIP_3D_Y" )
        write_3D_DATA( 0, state_k, EDIP_0k[:,:,:,2], "EL_DIP_3D_Z" )

        # Write MAGNETIC dipole density to cube file
        write_3D_DATA( 0, state_k, MDIP_0k[:,:,:,0], "MAG_DIP_3D_X" )
        write_3D_DATA( 0, state_k, MDIP_0k[:,:,:,1], "MAG_DIP_3D_Y" )
        write_3D_DATA( 0, state_k, MDIP_0k[:,:,:,2], "MAG_DIP_3D_Z" )

        # ELECTRIC dipole density (r,r')
        EDIP_0k_X = np.einsum( "xyzD->xD", EDIP_0k )
        EDIP_0k_Y = np.einsum( "xyzD->yD", EDIP_0k )
        EDIP_0k_Z = np.einsum( "xyzD->zD", EDIP_0k )
        write_2D_DATA( 0, state_k, np.einsum( "aD,xD->ax", EDIP_0k_X, EDIP_0k_X ), "EL_DIP_2D_X", 0 )
        write_2D_DATA( 0, state_k, np.einsum( "bD,yD->by", EDIP_0k_Y, EDIP_0k_Y ), "EL_DIP_2D_Y", 1 )
        write_2D_DATA( 0, state_k, np.einsum( "cD,zD->cz", EDIP_0k_Z, EDIP_0k_Z ), "EL_DIP_2D_Z", 2 )

        # MAGNETIC dipole density (r,r')
        MDIP_0k_X = np.einsum( "xyzD->xD", MDIP_0k )
        MDIP_0k_Y = np.einsum( "xyzD->yD", MDIP_0k )
        MDIP_0k_Z = np.einsum( "xyzD->zD", MDIP_0k )
        write_2D_DATA( 0, state_k, np.einsum( "aD,xD->ax", MDIP_0k_X, MDIP_0k_X ), "MAG_DIP_2D_X", 0 )
        write_2D_DATA( 0, state_k, np.einsum( "bD,yD->by", MDIP_0k_Y, MDIP_0k_Y ), "MAG_DIP_2D_Y", 1 )
        write_2D_DATA( 0, state_k, np.einsum( "cD,zD->cz", MDIP_0k_Z, MDIP_0k_Z ), "MAG_DIP_2D_Z", 2 )




        # Compute CD Density
        """
        \\vec{MU}_EL  = \int~dx dy dz  \\vec{x, y, z} \\rho_{TD}(x, y, z)
        \\vec{MU}_MAG = \int~dx'dy'dz' \\vec{x',y',z'} x \\nabla \\rho_{TD}(x',y',z')
        CD_0k         = \\vec{MU}_EL \cdot \\vec{MU}_MAG

        \\vec{MU}_EL(x, y, z)   = \\vec{x, y, z} \\rho_{TD}(x, y, z)
        \\vec{MU}_MAG(x',y',z') = \\vec{x',y',z'} x \\nabla \\rho_{TD}(x',y',z')
        CD_0k_2D(r,r') = \\vec{MU}_EL(x, y, z) * \\vec{MU}_MAG(x',y',z')
        CD_0k_EL(r)    = \int~dx'dy'dz' \\vec{MU}_EL(x, y, z) * \\vec{MU}_MAG(x',y',z')
        CD_0k_MAG(r')  = \int~dx dy dz  \\vec{MU}_EL(x, y, z) * \\vec{MU}_MAG(x',y',z')
        """


        # Total CD Density (r,r')
        EDIP_0k_X = np.einsum( "xyzD->xD", EDIP_0k )
        EDIP_0k_Y = np.einsum( "xyzD->yD", EDIP_0k )
        EDIP_0k_Z = np.einsum( "xyzD->zD", EDIP_0k )
        MDIP_0k_X = np.einsum( "xyzD->xD", MDIP_0k )
        MDIP_0k_Y = np.einsum( "xyzD->yD", MDIP_0k )
        MDIP_0k_Z = np.einsum( "xyzD->zD", MDIP_0k )
        write_2D_DATA( 0, state_k, np.einsum( "aD,xD->ax", EDIP_0k_X, MDIP_0k_X ), "CD_2D_X", 0 )
        write_2D_DATA( 0, state_k, np.einsum( "bD,yD->by", EDIP_0k_Y, MDIP_0k_Y ), "CD_2D_Y", 1 )
        write_2D_DATA( 0, state_k, np.einsum( "cD,zD->cz", EDIP_0k_Z, MDIP_0k_Z ), "CD_2D_Z", 2 )
        

        # Total Electric CD Density (r)
        MDIP_0k_T = np.einsum( "xyzD->D", MDIP_0k )
        write_3D_DATA( 0, state_k, np.einsum( "xyzD,D->xyz", EDIP_0k, MDIP_0k_T ), "CD_EL_3D_X" )
        write_3D_DATA( 0, state_k, np.einsum( "xyzD,D->xyz", EDIP_0k, MDIP_0k_T ), "CD_EL_3D_Y" )
        write_3D_DATA( 0, state_k, np.einsum( "xyzD,D->xyz", EDIP_0k, MDIP_0k_T ), "CD_EL_3D_Z" )

        # Total Electric CD Density (r')
        EDIP_0k_T = np.einsum( "xyzD->D", EDIP_0k )
        write_3D_DATA( 0, state_k, np.einsum( "D,xyzD->xyz", EDIP_0k_T, MDIP_0k ), "CD_MAG_3D_X" )
        write_3D_DATA( 0, state_k, np.einsum( "D,xyzD->xyz", EDIP_0k_T, MDIP_0k ), "CD_MAG_3D_Y" )
        write_3D_DATA( 0, state_k, np.einsum( "D,xyzD->xyz", EDIP_0k_T, MDIP_0k ), "CD_MAG_3D_Z" )








def main():
    get_Globals()
    compute_CD_Density()


if ( __name__ == "__main__" ):
    main()