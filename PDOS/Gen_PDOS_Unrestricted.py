import numpy as np
from matplotlib import pyplot as plt
import subprocess as sp

NFrags = 3
Npts = 2000
Emin = -5 # eV
Emax =  5 # eV
sig = 0.05 # eV

frag_list  = ["Se","Cd","Cl"]
color_list = ["black", "red", "green", "blue"] # FIRST COLOR IS TOTAL 


def getHOMOLUMOGap():
    lines = open('geometry.out','r').readlines()
    sp.call(' grep "Alpha  occ." geometry.out | tail -n 2 > occ_ALPHA_HOMO.dat',shell=True)
    sp.call(' grep "Gap" geometry.out | tail -n 2 | head -n 1 > ALPHA_HL_Gap.dat',shell=True)
    sp.call(' grep "Beta  occ." geometry.out | tail -n 1 > occ_BETA_HOMO.dat',shell=True)
    sp.call(' grep "Gap" geometry.out | tail -n 1 > BETA_HL_Gap.dat',shell=True)
    EHOMO_ALPHA = float( open('occ_ALPHA_HOMO.dat','r').readlines()[0].split()[-1] ) * 27.2114
    EHOMO_BETA  = float( open('occ_BETA_HOMO.dat','r').readlines()[0].split()[-1] ) * 27.2114
    EHOMO_TOTAL = np.max( [EHOMO_ALPHA, EHOMO_BETA]  )
    HLGap_ALPHA = float( open('ALPHA_HL_Gap.dat','r').readlines()[0].split()[1] ) * 27.2114
    HLGap_BETA = float( open('ALPHA_HL_Gap.dat','r').readlines()[0].split()[1] ) * 27.2114
    ELUMO_ALPHA = HLGap_ALPHA + EHOMO_ALPHA
    ELUMO_BETA = HLGap_BETA + EHOMO_BETA
    ELUMO_TOTAL = np.min( [ELUMO_ALPHA, ELUMO_BETA]  )
    E_Fermi_ALPHA = (ELUMO_ALPHA + EHOMO_ALPHA) / 2
    E_Fermi_BETA  = (ELUMO_BETA + EHOMO_BETA) / 2
    E_Fermi_TOTAL = (ELUMO_TOTAL + EHOMO_TOTAL) / 2
    HLGap_TOTAL = ELUMO_TOTAL - EHOMO_TOTAL
    print (f'\tE(HOMO) (ALPHA, BETA) = \t({round(EHOMO_ALPHA,4)}, {round(EHOMO_BETA,4)}) eV')
    print (f'\tE(HL)   (ALPHA, BETA, TOTAL) = \t({round(HLGap_ALPHA,4)}, {round(HLGap_BETA,4)}, {round(HLGap_TOTAL,4)}) eV')
    print (f'\tE(FERMI) (ALPHA, BETA) = \t({round(E_Fermi_ALPHA,4)}, {round(E_Fermi_BETA,4)}, {round(E_Fermi_TOTAL,4)}) eV')
    return E_Fermi_ALPHA, E_Fermi_BETA, E_Fermi_TOTAL

E_Fermi_ALPHA, E_Fermi_BETA, E_Fermi_TOTAL = getHOMOLUMOGap()


### PLOT SUMMED ALPHA AND BETA CONTRIBUTIONS ###

data = open('DOS_line_summSPINS.dat','r').readlines()[10:]
data = np.array([ np.array(j.split()).astype(float) for j in data if j.split()[0] != '**********' ])
Ei = data[:,0]
TDOS = data[:,1]
PDOS = data[:,3:3+NFrags]

PDOS_broaden = np.zeros(( Npts,NFrags+1 )) # INCLUDE TOTAL DOS IN FIRST COLUMN

dE = (Emax - Emin)/Npts
for n in range( Npts ):
    x = Emin + n*dE
    for frag in range( NFrags ):
        PDOS_broaden[n,1+frag] += np.sum( PDOS[:,frag] * np.exp( -( x - (Ei[:] - E_Fermi_TOTAL) )**2 / 2 / sig**2 ) )

norm = np.max( PDOS_broaden )
PDOS_broaden /= norm
for frag in range( NFrags ):
    plt.plot( np.arange(Emin,Emax,dE), PDOS_broaden[:,frag+1], c=color_list[frag+1], linewidth=2, label=frag_list[frag] )
plt.plot( np.arange(Emin,Emax,dE), np.sum(PDOS_broaden[:,:],axis=-1), color_list[0], linewidth=3, alpha=0.6, label="Total" )

plt.legend()
plt.xlim(Emin,Emax)
plt.ylim(0)
plt.xlabel('MO Energy (eV)',fontsize=15)
plt.ylabel('Projected DOS (Arb. Units)',fontsize=15)
#plt.title('Projected DOS',fontsize=15)
plt.savefig('PDOS_sumSPINS.jpg')
plt.clf()

PDOS_broaden[:,0] = np.sum( PDOS_broaden[:,:], axis=-1 )
np.savetxt('PDOS_sumSPINS.dat', np.c_[ np.arange(Emin,Emax,dE), PDOS_broaden ])




### PLOT ALPHA AND BETA CONTRIBUTIONS AS POSITIVE AND NEGATIVE ###

data_ALPHA = open('DOS_line_ALPHA.dat','r').readlines()[10:]
data_ALPHA = np.array([ np.array(j.split()).astype(float) for j in data_ALPHA if j.split()[0] != '**********' ])
Ei_ALPHA = data_ALPHA[:,0]
TDOS_ALPHA = data_ALPHA[:,1]
PDOS_ALPHA = data_ALPHA[:,3:3+NFrags]

PDOS_broaden_ALPHA = np.zeros(( Npts,NFrags+1 )) # INCLUDE TOTAL IN FIRST COLUMN


data_BETA = open('DOS_line_BETA.dat','r').readlines()[10:]
data_BETA = np.array([ np.array(j.split()).astype(float) for j in data_BETA if j.split()[0] != '**********' ])
Ei_BETA = data_BETA[:,0]
TDOS_BETA = data_BETA[:,1]
PDOS_BETA = data_BETA[:,3:3+NFrags]

PDOS_broaden_BETA = np.zeros(( Npts,NFrags+1 )) # INCLUDE TOTAL IN FIRST COLUMN

dE = (Emax - Emin)/Npts
for n in range( Npts ):
    x = Emin + n*dE
    for frag in range( NFrags ):
        PDOS_broaden_ALPHA[n,frag+1] += np.sum( PDOS_ALPHA[:,frag] * np.exp( -( x - (Ei_ALPHA[:] - E_Fermi_TOTAL) )**2 / 2 / sig**2 ) )
        PDOS_broaden_BETA[n,frag+1] += np.sum( PDOS_BETA[:,frag] * np.exp( -( x - (Ei_BETA[:] - E_Fermi_TOTAL) )**2 / 2 / sig**2 ) )

norm_ALPHA = np.max(PDOS_broaden_ALPHA)
PDOS_broaden_ALPHA /= norm_ALPHA
PDOS_broaden_BETA /= norm_ALPHA # NORMALIZE BOTH BY SAME REFERENCE

for frag in range( NFrags ):
    plt.plot( np.arange(Emin,Emax,dE),      PDOS_broaden_ALPHA[:,frag+1], c=color_list[frag+1], linewidth=3, label=frag_list[frag] )
    plt.plot( np.arange(Emin,Emax,dE), -1 * PDOS_broaden_BETA[:,frag+1],  c=color_list[frag+1], linewidth=3 )#, label=frag_list[frag] )

plt.plot( np.arange(Emin,Emax,dE), np.sum(PDOS_broaden_ALPHA[:,:],axis=-1), c=color_list[0], linewidth=3, alpha=0.5, label="Total" )
plt.plot( np.arange(Emin,Emax,dE), -1 * np.sum(PDOS_broaden_BETA[:,:],axis=-1),  c=color_list[0], linewidth=3, alpha=0.5 )

plt.legend()
plt.xlim(Emin,Emax)
#plt.ylim(0)
plt.xlabel('MO Energy (eV)',fontsize=15)
plt.ylabel('Projected DOS (Arb. Units)',fontsize=15)
#plt.title('Projected DOS',fontsize=15)
plt.savefig('PDOS_posALPHA_negBETA.jpg')
plt.clf()

PDOS_broaden_ALPHA[:,0] = np.sum( PDOS_broaden_ALPHA[:,:], axis=-1 )
PDOS_broaden_BETA[:,0] = np.sum( PDOS_broaden_BETA[:,:], axis=-1 )
np.savetxt('PDOS_ALPHA.dat', np.c_[ np.arange(Emin,Emax,dE), PDOS_broaden_ALPHA ])
np.savetxt('PDOS_BETA.dat', np.c_[ np.arange(Emin,Emax,dE), PDOS_broaden_BETA ])




### FOR FUTURE, SHOULD MAKE TOTAL DOS PLOT ON SAME PANEL FOR ALPHA AND BETA .... ###














