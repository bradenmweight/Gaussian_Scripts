import numpy as np
from matplotlib import pyplot as plt
import subprocess as sp


frag_list = ['C','Cl'] # This needs to be changed manually
Npts = 2000
Emin = -2 # eV
Emax = 2 # eV
EGRID = np.linspace(Emin,Emax,Npts)
sig = 0.05 # eV
STICKS   = True 
COLOR_IN = True



def getHOMOLUMOGap():
    lines = open('geometry.out','r').readlines()
    sp.call(' grep "Alpha  occ." geometry.out | tail -n 1 > occ_HOMO.dat',shell=True)
    sp.call(' grep "Gap" geometry.out | tail -n 1 > HL_Gap.dat',shell=True)
    EHOMOm1 = float( open('occ_HOMO.dat','r').readlines()[0].split()[-2] ) * 27.2114
    EHOMO = float( open('occ_HOMO.dat','r').readlines()[0].split()[-1] ) * 27.2114
    HLGap = float( open('HL_Gap.dat','r').readlines()[0].split()[1] ) * 27.2114
    ELUMO = HLGap + EHOMO
    EFERMI = (ELUMO + EHOMO)/2
    print (f'\tE(HOMO-1) = {EHOMOm1-EFERMI} eV')
    print (f'\tE(HOMO) = {EHOMO-EFERMI} eV')
    print (f'\tE(LUMO) = {ELUMO-EFERMI}')
    print (f'\tE(HL) = {HLGap} eV')
    print (f'\tE(FERMI) = {EFERMI} eV')
    return EFERMI, HLGap



E_Fermi, HLGap = getHOMOLUMOGap()
NFrags = len(frag_list)
data = open('DOS_line.txt','r').readlines()[36:]
data = np.array([ np.array(j.split()).astype(float) for j in data if j.split()[0] != '**********' ])
data = np.array([ j for j in data if not np.allclose(j[1:], np.zeros(len(j[1:]))) ])
Ei = data[:,0]
TDOS = data[:,1]
PDOS = data[:,3:3+NFrags] / np.max(TDOS)

PDOS_broaden = np.zeros(( Npts,NFrags+1 ))
TDOS_broaden = np.zeros(( Npts ))

dE = (Emax - Emin)/Npts
for n in range( Npts ):
    x = Emin + n*dE
    for frag in range( NFrags ):
        PDOS_broaden[n,frag] += np.sum( PDOS[:,frag] * np.exp( -( x - (Ei[:] - E_Fermi) )**2 / 2 / sig**2 ) )
    TDOS_broaden[n] += np.sum( 1.000 * np.exp( -( x - (Ei[:] - E_Fermi) )**2 / 2 / sig**2 ) )
PDOS_broaden[:,-1] = np.sum( PDOS_broaden[:,:-1],axis=-1 )

#norm = np.max( PDOS_broaden )
#PDOS_broaden /= norm
color_list = ['blue', 'red', 'black', 'green', 'orange']
for frag in range( NFrags ):
    plt.plot( np.arange(Emin,Emax,dE), PDOS_broaden[:,frag], c=color_list[frag], label=f'{frag_list[frag]}')
    if (STICKS == True) :
        plt.scatter( Ei - E_Fermi, PDOS[:,frag], s=30, c=color_list[frag], alpha=0.5, edgecolors=color_list[frag] )
    if ( COLOR_IN == True ):
        plt.fill_between( x=EGRID, y1=PDOS_broaden[:,frag], where = (EGRID < 0), color=color_list[frag], alpha=0.1 )
#plt.plot( np.arange(Emin,Emax,dE), PDOS_broaden[:,-1], c='black', label=f'Total')


plt.legend(loc='upper center')
plt.xlim(Emin,Emax)
plt.ylim(-0.01)
plt.title('$E_{HL}$'+' = %1.2f eV'%(HLGap),fontsize=15)
plt.xlabel('MO Energy (eV)',fontsize=15)
plt.ylabel('Projected DOS (Arb. Units)',fontsize=15)
plt.savefig('PDOS.jpg', dpi=300)

np.savetxt('PDOS.dat', np.c_[ np.linspace(Emin,Emax,Npts), PDOS_broaden ])

