import numpy as np
import matplotlib.pyplot as plt

from pyscf import gto, dft, tddft
from pyscf.tools import cubegen

def get_atomic_numbers( mol ):
    return [ mol.atom_charge(i) for i in range(mol.natm) ]

def get_electric_dipole_ao( mol, Shalf=None ):
    charges    = mol.atom_charges()
    coords     = mol.atom_coords()
    nuc_dipole = np.einsum("a,ad->d", charges, coords) / charges.sum()
    with mol.with_common_orig(nuc_dipole):
        dipole_ao  = mol.intor_symmetric("int1e_r", comp=3)
    if ( Shalf is not None ):
        dipole_ao = np.einsum( 'ap,xpq,bq->xab', Shalf, dipole_ao, Shalf )
    return dipole_ao

def get_magnetic_dipole_ao( mol, Shalf=None ):
    charges    = mol.atom_charges()
    coords     = mol.atom_coords()
    nuc_dipole = np.einsum("a,ad->d", charges, coords) / charges.sum()
    with mol.with_common_orig(nuc_dipole):
        #dipole_ao  = -mol.intor_symmetric("int1e_cg_irxp", comp=3) # SYMMETRIC DOES NOT WORK FOR r x p
        dipole_ao  = -mol.intor("int1e_cg_irxp", comp=3)
    if ( Shalf is not None ):
        dipole_ao = np.einsum( 'ap,xpq,bq->xab', Shalf, dipole_ao, Shalf )
    return dipole_ao

def write_CUBE( ORIGIN, atom_labels, atom_numbers, coords, Lxyz, dLxyz, DATA, FILE_NAME ):
    NAtoms = len(coords)
    Nxyz = DATA.shape[:3]
    f = open(f'{FILE_NAME}.cube','w')
    f.write(f"Circular Dichroism Density \n")
    f.write(f"Totally {np.prod(Nxyz)} grid points\n")
    f.write(f"{NAtoms} {Lxyz[0]} {Lxyz[1]} {Lxyz[2]}\n")
    f.write(f"{Nxyz[0]} {dLxyz[0]}  0.000000   0.000000\n")
    f.write(f"{Nxyz[1]} 0.000000   {dLxyz[1]} 0.000000\n")
    f.write(f"{Nxyz[2]} 0.000000   0.000000   {dLxyz[2]} \n")
    for at in range(NAtoms):
        #f.write( "%d %1.6f %1.6f %1.6f %1.6f\n" % (int(atom_numbers[at]),atom_numbers[at],coords[at,0]-ORIGIN[0],coords[at,1]-ORIGIN[1],coords[at,2]-ORIGIN[2]) )
        f.write( "%d %1.6f %1.6f %1.6f %1.6f\n" % (int(atom_numbers[at]),atom_numbers[at],coords[at,0],coords[at,1],coords[at,2]) )
    for x in range(Nxyz[0]):
        #print(f'X = {x}')
        for y in range(Nxyz[1]):
            outArray = []
            for z in range(Nxyz[2]):
                outArray.append( DATA[x,y,z] )
                if ( len(outArray) % 6 == 0 or z == Nxyz[2]-1 ):
                    #outArray.append('\n')
                    f.write( " ".join(map( str, np.round(outArray,8) )) + "\n" )
                    outArray = []
    f.close()

def main():
    mol = gto.Mole()
    mol.atom = """
 C                  0.72599764   -1.17622359   -0.00284968
 C                  2.12115764   -1.17622359   -0.00284968
 C                  2.81869564    0.03152741   -0.00284968
 C                  2.12104164    1.24003641   -0.00404868
 C                  0.72621664    1.23995841   -0.00452768
 C                  0.02861564    0.03175241   -0.00353168
 H                  0.17623864   -2.12854059   -0.00239968
 H                  2.67066564   -2.12873659   -0.00153468
 H                  2.67124164    2.19217941   -0.00410768
 H                  0.17609464    2.19223941   -0.00548068
 H                 -1.07098836    0.03193541   -0.00371168
 C                  4.35869538    0.03163944   -0.00196182
 H                  4.71594233    0.04405522   -1.01048573
 F                  4.80824061    1.12604337    0.64818366
 Cl                 4.94499941   -1.41543325    0.81033628
    """
    mol.basis = "sto-3g"
    mol.unit = "angstrom"
    mol.verbose = 0
    mol.symmetry = False
    mol.build()

    mf = dft.RKS(mol)
    mf.xc = "b3lyp"
    mf.kernel()

    td = tddft.TDDFT(mf)
    td.nstates = 1
    td.kernel()

    print("Excitation energies (eV):")
    print(td.e * 27.2114)
    print("Oscillator strengths:")
    print(td.oscillator_strength())
    # Compute the electric dipole transition moments
    electric_dipole = td.transition_dipole()
    print("PYSCF Electric dipole transition moments (a.u.):")
    print(electric_dipole)
    # Compute the magnetic dipole transition moments
    magnetic_dipole = td.transition_magnetic_dipole()
    print("PYSCF Magnetic dipole transition moments (a.u.):")
    print(magnetic_dipole)

    n_elec_alpha, n_elec_beta = mol.nelec
    S_AO       = mol.intor('int1e_ovlp') # AO overlaps
    N_AO       = len(S_AO)
    occ_inds, vir_inds = slice( n_elec_alpha ), slice( n_elec_alpha, N_AO )
    C_MO       = mf.mo_coeff # MO coefficients in non-orthogonal AO basis
    XY         = np.array(td.xy) # Excitation/de-excitation TDDFT coefficients (NSTATES,(x,y),OCC,VIR)
    X,Y        = XY[:,0], XY[:,1] # Excitation/de-excitation TDDFT coefficients (NSTATES,OCC,VIR)
    el_dip_ao  = get_electric_dipole_ao( mol )
    el_dip_MO = np.einsum( 'ai,Dab,bj->Dij', C_MO.conj(), el_dip_ao, C_MO ) # Dipole matrix in MO basis
    el_dip    = 2 * np.einsum( 'Jij,Dij->JD', X+Y, el_dip_MO[:,occ_inds,vir_inds] ) # Dipole matrix in ES
    print("BRADEN Electric dipole transition moments (a.u.):")
    print(el_dip)
    mag_dip_ao = get_magnetic_dipole_ao( mol )
    mag_dip_MO = np.einsum( 'ai,Dab,bj->Dij', C_MO.conj(), mag_dip_ao, C_MO ) # Dipole matrix in MO basis
    mag_dip    = 2 * np.einsum( 'Jij,Dij->JD', X-Y, mag_dip_MO[:,occ_inds,vir_inds] ) # Dipole matrix in ES
    print("BRADEN Magnetic dipole transition moments (a.u.):")
    print(mag_dip)


    # Project the transition density X into real-space grid and sace as a cube file
    nx,ny,nz = 50,50,50
    ORIGIN   = np.average( mol.atom_coords(), axis=0 )
    MIN_X    = np.min( mol.atom_coords()[:,0] ) - 1
    MAX_X    = np.max( mol.atom_coords()[:,0] ) + 1
    MIN_Y    = np.min( mol.atom_coords()[:,1] ) - 1
    MAX_Y    = np.max( mol.atom_coords()[:,1] ) + 1
    MIN_Z    = np.min( mol.atom_coords()[:,2] ) - 1
    MAX_Z    = np.max( mol.atom_coords()[:,2] ) + 1
    Lxyz     = np.array( [MAX_X-MIN_X, MAX_Y-MIN_Y, MAX_Z-MIN_Z] )
    dLxyz    = Lxyz / np.array( [nx,ny,nz] )
    grid     = np.meshgrid( np.linspace(MIN_X,MAX_X,nx), 
                            np.linspace(MIN_Y,MAX_Y,ny), 
                            np.linspace(MIN_Z,MAX_Z,nz) )
    grid     = np.array( grid ).reshape( (3,-1) ).T # {|x,y,z>}
    AO_xyz   = mol.eval_gto("GTOval_cart", coords=grid  ).reshape( (nx,ny,nz,-1) ) # < x,y,z | AO_i >

    # Transform transition density X into AO and then to real-space grid
    X_AO     = np.einsum( 'ai,Jij,bj->Jab', C_MO[:,occ_inds].conj(), X, C_MO[:,vir_inds] ) # Transition density in AO basis
    Y_AO     = np.einsum( 'ai,Jij,bj->Jab', C_MO[:,occ_inds].conj(), Y, C_MO[:,vir_inds] ) # Transition density in AO basis
    X_xyz    = np.einsum( 'xyza,Jab,xyzb->Jxyz', AO_xyz, X_AO, AO_xyz )
    Y_xyz    = np.einsum( 'xyza,Jab,xyzb->Jxyz', AO_xyz, Y_AO, AO_xyz )
    XpY_xyz  = X_xyz + Y_xyz


    HOMO_xyz = np.einsum( 'a,xyza->xyz', C_MO[:,n_elec_alpha-1], AO_xyz )

    atom_labels  = [ mol.atom_symbol(i) for i in range(mol.natm) ]
    atom_numbers = get_atomic_numbers( mol )
    #for state in range( X.shape[0] ):
        #myCUBE = cubegen.Cube( mol, nx=nx, ny=ny, nz=nz )
        #myCUBE = cubegen.Cube( mol, nx=nx, ny=ny, nz=nz, margin=0.0, origin=ORIGIN, extent=[MIN_X,MAX_X,MIN_Y,MAX_Y,MIN_Z,MAX_Z] )
        #myCUBE.write( field=XpY_xyz[state], fname="TD_0_%d.cube" % (state+1) )
        #write_CUBE( ORIGIN, atom_labels, atom_numbers, mol.atom_coords(), Lxyz, dLxyz, XpY_xyz[state], f"TD_0_{state+1}" )
    write_CUBE( ORIGIN, atom_labels, atom_numbers, mol.atom_coords(), Lxyz, dLxyz, HOMO_xyz, f"HOMO" )




if ( __name__ == "__main__" ):
    main()
