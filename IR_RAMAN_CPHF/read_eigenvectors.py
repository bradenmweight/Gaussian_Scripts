import numpy as np
import os
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy.interpolate import interp1d

COM_FILE = open("geometry.com","r").readlines()
FCHK_FILE = open("geometry.fchk","r").readlines()

NAtoms = 418

at_labels = []
geom      = np.zeros(( NAtoms, 3 ))
for count, line in enumerate( COM_FILE ):
    if ( "0 1".split() == line.split() ):
        for at in range( NAtoms ):
            t = COM_FILE[ count+1+at ].split()
            at_labels.append( t[0] )
            geom[at,:] = np.array( t[1:] ).astype(float)

if ( not os.path.isfile("NM_EIGV.dat.npy") ):
    NM_EIGV = []
    for count, line in enumerate( FCHK_FILE ):
        t = line.split()
        if ( "Vib-Modes" in t ):
            length = int( t[-1] )
            NColumns = len(FCHK_FILE[ count+1 ].split())
            for n in range( int(np.ceil(length/NColumns)) ):
                s = FCHK_FILE[ count+1+n ].split()
                for j in range( len(s) ):
                    NM_EIGV.append( float(s[j]) )

    NM_EIGV = np.array( NM_EIGV )
    print(len(NM_EIGV))
    NModes = len(NM_EIGV) // NAtoms // 3
    print(NModes)
    NM_EIGV = NM_EIGV.reshape(( NModes, NAtoms, 3 ))
    np.save( "NM_EIGV.dat", NM_EIGV )
else:
    NM_EIGV = np.load( "NM_EIGV.dat.npy" )
    NModes = np.shape(NM_EIGV)[0]

# Check normal mode modes between off- and on-resonant cases
RES_NM = np.load("../CPHF_SCAN/1/NM_EIGV.npy")
count = 0
for m in range( len(NM_EIGV) ):
    for at in range( NAtoms ):
        if ( (np.round(NM_EIGV[m,at,:],2) != RES_NM[m,at,:]).all() ):
            print( "Off-resonant: atom, vec =", at, np.round(NM_EIGV[m,at,:],2) )
            print( " On-resonant: atom, vec =", at, RES_NM[m,at,:] )
            count += 1
    if( count > 0 ):
        print( f"Mode Index {m} did not match" )
        break
if( count == 0 ):
    print( "Eigenvectors match!" )
else:
    print( "Eigenvectors do not match!" )
    exit()





# Find SWCNT axis
L = [0,0,0]
L[0] = np.max( geom[:,0] ) - np.min( geom[:,0] )
L[1] = np.max( geom[:,1] ) - np.min( geom[:,1] )
L[2] = np.max( geom[:,2] ) - np.min( geom[:,2] )
print(f"Tube is {L[np.argmax(L)]} "+f"$\AA$"+" along axis {np.argmax(L)}")
Lmain = L[np.argmax(L)]
dmain = np.argmax(L)
L_min = np.min( geom[:,dmain] )
L_max = np.max( geom[:,dmain] )

geom[:,dmain] = geom[:,dmain] - L_min

# Bin vector magnitudes along SWCNT axis
NBins = 100
NM_EIGV_norm = np.linalg.norm( NM_EIGV[:,:,:], axis=-1 )

BINS = np.linspace(0, Lmain, NBins )
dBin = Lmain / NBins 
mode_localization_hist = np.zeros(( NModes, NBins ))

if ( not os.path.isfile(f"mode_localization_hist_NBins{NBins}.dat.npy") ):
    # This is most expensive part...Probably could use numpy histogram function
    for b in range( NBins-1 ):
        print(b)
        for m in range( NModes ):
            for at in range( NAtoms ):
                if ( BINS[b] <= geom[at,dmain] and geom[at,dmain] < BINS[b+1] ):
                    mode_localization_hist[m,b] += NM_EIGV_norm[m,at]
    np.save(f"mode_localization_hist_NBins{NBins}.dat", mode_localization_hist)
else:
    mode_localization_hist = np.load( f"mode_localization_hist_NBins{NBins}.dat.npy" )


# Normalize probability for each mode
for m in range( NModes ):
    mode_localization_hist[m,:] /= np.sum( mode_localization_hist[m,:] )

# Read frequencies to convolve with populations. This is generated from ~/Gaussain_scripts/IR_RAMAN_CPHF/get_IR_Raman_SPEC.py
FREQ_IR_RAMAN = np.loadtxt("FREQ_IR_RAMAN.dat")
FREQ = FREQ_IR_RAMAN[:,0]
cmap = mpl.cm.magma

# Plot heatmap as function of normal mode index
plt.imshow( np.sqrt(mode_localization_hist[:,:]), cmap=cmap , origin="lower", aspect='auto', extent=[BINS[0], BINS[-1], 0, NModes-1], vmin=0, vmax=0.8 )
plt.xlabel(r"SWCNT Axis Position ($\AA$)",fontsize=15)
plt.ylabel("Normal Mode Index",fontsize=15)
plt.colorbar( pad=0.01 )
plt.tight_layout()
plt.savefig(f"Mode_HEATMAP_NBINS{NBins}.jpg", dpi=500)
plt.clf()
np.savetxt( f"Mode_HEATMAP_NBINS{NBins}_AMPLITUDE.dat", np.sqrt(mode_localization_hist[:,:]) )

# Plot heatmap as function of normal mode frequency
NPTS = 2000
sig = 5 # cm^-1
EMIN = 0.0  # cm^-1
EMAX = 2000 # cm^-1
dE = ( EMAX - EMIN ) / NPTS
SPEC = np.zeros(( NPTS, NBins ))
SPEC_mode = np.zeros(( NModes, NBins ))
for pt in range( NPTS ):
    E = EMIN + dE * pt
    for b in range( NBins ):
        SPEC[pt,b] = np.sum( np.sqrt(mode_localization_hist[:,b]) * np.exp( -(E - FREQ[:])**2 / 2 / sig**2 ) )
        #SPEC_mode[:,b] = mode_localization_hist[:,b] * np.exp( -(E - FREQ[:])**2 / 2 / (sig/100)**2 )

plt.imshow( SPEC, cmap=cmap , origin="lower", aspect='auto', extent=[BINS[0], BINS[-1], FREQ[0], FREQ[-1]])
plt.xlim( 0, BINS[-1] )
plt.ylim( 0, EMAX )
#plt.ylim( 0, 1000 )
#plt.ylim( 1000, 2000 )
plt.xlabel(r"SWCNT Axis Position ($\AA$)",fontsize=15)
plt.ylabel(r"Normal Mode Frequency (cm$^{-1}$)",fontsize=15)
plt.colorbar( pad=0.01 )
plt.tight_layout()
plt.savefig(f"Mode_HEATMAP_ENERGY_NBINS{NBins}.jpg", dpi=500)
plt.clf()

np.savetxt(f"Mode_HEATMAP_ENERGY_NBINS{NBins}.dat", SPEC)



# Plot some histograms

for m in range( NModes ):
    plt.plot( BINS, mode_localization_hist[m,:], label=f"Mode {m+1}" )
#plt.legend()
plt.xlabel("SWCNT Axis Position",fontsize=15)
plt.ylabel(r"Probability, P$_m$(X)",fontsize=15)
plt.tight_layout()
plt.savefig(f"Mode_Hist_NBINS{NBins}.jpg", dpi=500)
plt.clf()


"""
# Plot some histograms decomposed

for m in range( 1000 ):
    plt.plot( BINS, SPEC_mode[m,:], label=f"Mode {m+1}" )
#plt.legend()
plt.xlabel("SWCNT Axis Position",fontsize=15)
plt.ylabel(r"Probability, P$_m$(X)",fontsize=15)
plt.tight_layout()
plt.savefig(f"Mode_Hist_decomposed_NBINS{NBins}.jpg", dpi=500)
plt.clf()
"""










# ANALYZE EDGE STATES WITH OFF-RESONANT SPECTRA
MAX_IR    = np.max( FREQ_IR_RAMAN[:,1] )
MAX_RAMAN = np.max( FREQ_IR_RAMAN[:,2] )
N_AROUND_EDGE = int(5 / dBin)//2 # Let's count outer +-0.5 nm = +-5 \AA
print("Edges Index, N_Around_Edge, dBin =", [0,-1], N_AROUND_EDGE, dBin )
EDGE_DATA = []
for m in range( NModes ):
    if ( np.sum(mode_localization_hist[m,:N_AROUND_EDGE]) + np.sum(mode_localization_hist[m,-N_AROUND_EDGE:]) > 0.35 ): # If localized to edge more than 50 %
        EDGE_DATA.append( [m, FREQ_IR_RAMAN[m,0], FREQ_IR_RAMAN[m,1], FREQ_IR_RAMAN[m,2], np.sum(mode_localization_hist[m,:N_AROUND_EDGE]) + np.sum(mode_localization_hist[m,-N_AROUND_EDGE:])*100, FREQ_IR_RAMAN[m,1]/MAX_IR * 100, FREQ_IR_RAMAN[m,2]/MAX_RAMAN*100] )

EDGE_DATA = np.array(EDGE_DATA)
np.savetxt( f"EDGE_DATA_NBINS{NBins}.dat", np.c_[EDGE_DATA[:,0], EDGE_DATA[:,1], EDGE_DATA[:,2], EDGE_DATA[:,3], EDGE_DATA[:,4], EDGE_DATA[:,5], EDGE_DATA[:,6] ], header="Index, FREQ, IR, RAMAN, %EDGE, %IR_max, %RAMAN_max" )

# Plot edge data against frequency as percent of IR and RAMAN for system
plt.plot( EDGE_DATA[:,0], EDGE_DATA[:,-2], "o", label="IR / MAX(IR)" )
plt.plot( EDGE_DATA[:,0], EDGE_DATA[:,-1], "o", label="RAMAN / MAX(RAMAN)" )
plt.legend()
plt.xlim( 0, 2000 )
plt.ylim( 0)
plt.xlabel(r"Normal Mode Frequency (cm$^{-1}$)",fontsize=15)
plt.ylabel(r"Percent of Max Off-resonant Intensity",fontsize=15)
plt.title("Edge Contribution to off-resonant IR and Raman",fontsize=15)
plt.tight_layout()
plt.savefig(f"Edge_Mode_Percent_NBINS{NBins}.jpg", dpi=500)
plt.clf()



# ANALYZE EDGE STATES WITH OFF-RESONANT SPECTRA
CPHF_JOBS = ["1", "2", "3", "4", "5", "6", "7"]
CPHF_JOBS_FREQ = ["1.05 eV", "1.15 eV", "1.20 eV", "1.30 eV", "1.35 eV", "1.45 eV", "1.75 eV"]
for ind, job in enumerate(CPHF_JOBS):
    FREQ      = np.loadtxt(f"../CPHF_SCAN/{job}/IR_RAMAN_RAW.dat")[:,0]
    RAMAN     = np.loadtxt(f"../CPHF_SCAN/{job}/IR_RAMAN_RAW.dat")[:,3]
    MAX_RAMAN = np.max( RAMAN )
    EDGE_DATA = []
    for m in range( NModes ):
        if ( FREQ_IR_RAMAN[m,0] > 2000 ):
            continue
        if ( np.sum(mode_localization_hist[m,:N_AROUND_EDGE]) + np.sum(mode_localization_hist[m,-N_AROUND_EDGE:]) > 0.35 ): # If localized to edge more than 50 %
            EDGE_DATA.append( [m, FREQ_IR_RAMAN[m,0], FREQ_IR_RAMAN[m,1], FREQ_IR_RAMAN[m,2], np.sum(mode_localization_hist[m,:N_AROUND_EDGE]) + np.sum(mode_localization_hist[m,-N_AROUND_EDGE:])*100, RAMAN[m]/MAX_RAMAN*100])
    EDGE_DATA = np.array(EDGE_DATA)
    #print("# of edge modes =", len(EDGE_DATA))
    np.savetxt( f"EDGE_DATA_CPHFdir{job}_NBINS{NBins}.dat", np.c_[EDGE_DATA[:,0], EDGE_DATA[:,1], EDGE_DATA[:,2], EDGE_DATA[:,3], EDGE_DATA[:,4], EDGE_DATA[:,5] ], header="Index, FREQ, RAMAN, %EDGE, %RAMAN_max" )

    # Plot edge data against frequency as percent resonant RAMAN for system
    plt.plot( EDGE_DATA[:,1], EDGE_DATA[:,-1], "o", label="RAMAN / MAX(RAMAN)" )
    plt.xlim( 0, 2000 )
    plt.ylim( 0)
    plt.xlabel(r"Normal Mode Frequency (cm$^{-1}$)",fontsize=15)
    plt.ylabel(r"Percent of Max Intensity",fontsize=15)
    plt.title(r"Edge Contribution to E$_{CPHF}$ = "+f"{CPHF_JOBS_FREQ[ind]}",fontsize=15)
    plt.tight_layout()
    plt.savefig(f"Edge_Mode_Percent_CPHFdir{job}_NBINS{NBins}.jpg", dpi=500)
    plt.clf()
















# ANALYZE DEFECT STATES WITH OFF-RESONANT SPECTRA
MAX_IR    = np.max( FREQ_IR_RAMAN[:,1] )
MAX_RAMAN = np.max( FREQ_IR_RAMAN[:,2] )
MIDDLE_INDEX = NBins // 2 
N_AROUND_CENTER = int(10 / dBin)//2 # Let's count center +-0.5 nm = +-5 \AA
print("Middle Index, N_Around, dBin =", MIDDLE_INDEX, N_AROUND_CENTER, dBin )
INIDICES = [ j for j in range( MIDDLE_INDEX - N_AROUND_CENTER, MIDDLE_INDEX + N_AROUND_CENTER) ]

DEFECT_DATA = []
for m in range( NModes ):
    if ( np.sum(mode_localization_hist[m,INIDICES]) > 0.35 ): # If localized to center more than 50 %
        DEFECT_DATA.append( [m, FREQ_IR_RAMAN[m,0], FREQ_IR_RAMAN[m,1], FREQ_IR_RAMAN[m,2], np.sum(mode_localization_hist[m,INIDICES])*100, FREQ_IR_RAMAN[m,1]/MAX_IR * 100, FREQ_IR_RAMAN[m,2]/MAX_RAMAN*100] )

DEFECT_DATA = np.array(DEFECT_DATA)
np.savetxt( f"Defect_Mode_Percent_NBINS{NBins}.dat", np.c_[DEFECT_DATA[:,0], DEFECT_DATA[:,1], DEFECT_DATA[:,2], DEFECT_DATA[:,3], DEFECT_DATA[:,4], DEFECT_DATA[:,5], DEFECT_DATA[:,6] ], header="Index, FREQ, IR, RAMAN, %CENTER, %IR_max, %RAMAN_max" )

# Plot defect data against frequency as percent of IR and RAMAN for system
plt.plot( DEFECT_DATA[:,0], DEFECT_DATA[:,-2], "o", label="IR / MAX(IR)" )
plt.plot( DEFECT_DATA[:,0], DEFECT_DATA[:,-1], "o", label="RAMAN / MAX(RAMAN)" )
plt.legend()
plt.xlim( 0, 2000 )
plt.ylim( 0)
plt.xlabel(r"Normal Mode Frequency (cm$^{-1}$)",fontsize=15)
plt.ylabel(r"Percent of Max Off-resonant Intensity",fontsize=15)
plt.title("Defect Contribution to off-resonant IR and Raman",fontsize=15)
plt.tight_layout()
plt.savefig(f"Defect_Mode_Percent_NBINS{NBins}.jpg", dpi=500)
plt.clf()




#ANALYZE DEFECT STATES WITH RESONANT SPECTRA
CPHF_JOBS = ["1", "2", "3", "4", "5", "6", "7"]
CPHF_JOBS_FREQ = ["1.05 eV", "1.15 eV", "1.20 eV", "1.30 eV", "1.35 eV", "1.45 eV", "1.75 eV"]
for ind, job in enumerate(CPHF_JOBS):
    FREQ      = np.loadtxt(f"../CPHF_SCAN/{job}/IR_RAMAN_RAW.dat")[:,0]
    RAMAN     = np.loadtxt(f"../CPHF_SCAN/{job}/IR_RAMAN_RAW.dat")[:,3]
    MAX_RAMAN = np.max( RAMAN )
    RAMAN_INTERP = interp1d( FREQ, RAMAN, kind='cubic' )
    EDGE_DATA = []
    for m in range( NModes ):
        if ( FREQ_IR_RAMAN[m,0] > 2000 ):
            continue
        if ( np.sum(mode_localization_hist[m,INIDICES]) > 0.35 ): # If localized to center more than 50 %
            #print( INIDICES, np.sum(mode_localization_hist[m,INIDICES]) )
            EDGE_DATA.append( [m, FREQ_IR_RAMAN[m,0], FREQ_IR_RAMAN[m,1], FREQ_IR_RAMAN[m,2], np.sum(mode_localization_hist[m,INIDICES])*100, RAMAN_INTERP(FREQ_IR_RAMAN[m,0])/MAX_RAMAN*100] )
    EDGE_DATA = np.array(EDGE_DATA)
    #print("# of defect modes =", len(EDGE_DATA))
    np.savetxt( f"Defect_Mode_Percent_CPHFdir{job}_NBINS{NBins}.dat", np.c_[EDGE_DATA[:,0], EDGE_DATA[:,1], EDGE_DATA[:,2], EDGE_DATA[:,3], EDGE_DATA[:,4], EDGE_DATA[:,5] ], header="Index, FREQ, RAMAN, %CENTER, %RAMAN_max" )

    # Plot edge data against frequency as percent resonant RAMAN for system
    plt.plot( EDGE_DATA[:,1], EDGE_DATA[:,-1], "o", label="RAMAN / MAX(RAMAN)" )
    plt.xlim( 0, 2000 )
    plt.ylim( 0)
    plt.xlabel(r"Normal Mode Frequency (cm$^{-1}$)",fontsize=15)
    plt.ylabel(r"Percent of Max Intensity",fontsize=15)
    plt.title(r"Defect Contribution to E$_{CPHF}$ = "+f"{CPHF_JOBS_FREQ[ind]}",fontsize=15)
    plt.tight_layout()
    plt.savefig(f"Defect_Mode_Percent_CPHFdir{job}_NBINS{NBins}.jpg", dpi=500)
    plt.clf()






# Compute delocalization parameter, Ld
Ld = np.zeros((NModes))
P  = np.zeros((NModes, NBins))
for m in range( NModes ):
    tmp = np.sum( mode_localization_hist[m,:] ** 2 )
    Ld[m] = 1/tmp * (dBin / Lmain) * 100 # Write as percent of total length of SWCNT
np.savetxt(f"Ld_NBINS{NBins}.dat", np.c_[ FREQ_IR_RAMAN[:,0], Ld[:] ], header="FREQ, Ld")

plt.plot( FREQ_IR_RAMAN[:,0], Ld[:], "o" )
plt.xlabel(r"Normal Mode Frequency (cm$^{-1}$)",fontsize=15)
plt.ylabel(r"Normal Mode Delocalization, L$_d$ ($\AA$)",fontsize=15)
plt.xlim(0,2000)
plt.savefig( f"Ld_NBINS{NBins}.jpg" , dpi=500)


outARRAY = np.zeros(( len(CPHF_JOBS)+2, len(FREQ_IR_RAMAN[:,0]) ))
outARRAY[0,:] = FREQ_IR_RAMAN[:,0]
outARRAY[1,:] = Ld[:]
for ind, job in enumerate(CPHF_JOBS):
    RAMAN     = np.loadtxt(f"../CPHF_SCAN/{job}/IR_RAMAN_RAW.dat")[:,3]
    outARRAY[ind+2,:] = RAMAN[:] / np.max(RAMAN) * 100
np.savetxt(f"Ld_allDATA_NBINS{NBins}.dat", outARRAY.T, header=" ".join(["1.05 eV", "1.15 eV", "1.20 eV", "1.30 eV", "1.35 eV", "1.45 eV", "1.75 eV"]))
