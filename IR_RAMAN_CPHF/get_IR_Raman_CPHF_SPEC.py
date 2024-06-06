import numpy as np
from matplotlib import pyplot as plt

# Read Data
filename = "geometry.out"
NAtoms = 396
f = open(filename,"r").readlines() # Open File
FREQ = []
IR_INT = []
RAMAN = []
RAMAN_f1 = [] # Optical Frequency of Zero
RAMAN_f2 = [] # First non-zero optical frequency
NM_EIGV = np.zeros(( 3 * NAtoms, NAtoms, 3 ))
Optical_freq = []
counter = 0
good = 0
bad = 0
for count, line in enumerate(f):
    t = line.split()
    if ( len(t) >= 5 ):
        if ( t[0] == "Incident" and t[1] == "light" and t[2] == "(cm**-1):" ):
            for j in range( len(t)-3 ):
                Optical_freq.append( t[j+3] )

    if ( len(t) == 5 ):
        if ( t[0] == "Frequencies" and t[1] == "--" ):
            FREQ.append( float(t[2]) ) # cm^-1
            FREQ.append( float(t[3]) )
            FREQ.append( float(t[4]) )
    if ( len(t) >= 4 ):
        if ( t[0] == "IR" and t[1] == "Inten" ):
            IR_INT.append( float(t[3]) ) # KM/Mole
            IR_INT.append( float(t[4]) )
            IR_INT.append( float(t[5]) )
        if ( t[0] == "Raman" and t[1] == "Activ" ):
            RAMAN.append( float(t[3]) ) # A**4 / AMU
            RAMAN.append( float(t[4]) )
            RAMAN.append( float(t[5]) )        

        if ( t[0] == "RamAct" and t[1] == "Fr=" and t[2] == "1--" ):
            for j in range(3):
                RAMAN_f1.append( float(t[3+j]) ) # A**4 / AMU     

        # Read in Raman activity at first non-zero optical frequency

        if ( t[0] == "RamAct" and t[1] == "Fr=" and t[2].split('--')[0] == "2" ):
            #print (t)
            if ( len(t) == 5 ):
                for j in range( 5 ):
                    if( j == 2 ):
                        t[j] = t[j].split('--')[1]
                    else:
                        t[j] = t[j]
                    t.append
            counter += 3
            for j in range(3):
                #print(j, len(RAMAN_f2))
                try:
                    RAMAN_f2.append( float(t[-3+j]) ) # A**4 / AMU
                    print("\tOKAY")
                    good += 1
                    #counter += 1
                except ValueError:
                    #counter += 1
                    print ("Found ValueError. Setting this intensity to zero")
                    RAMAN_f2.append( -1.0 )
                    bad += 1
                    #print(f"\t\t{len(RAMAN_f2)}")

        if ( t == "Atom AN X Y Z X Y Z X Y Z".split() ):
            mode_IND = len(RAMAN_f2)-1
            for at in range( NAtoms ):
                NM_EIGV[mode_IND-2,at,:] = np.array( f[count+1 + at].split()[2:5] ).astype(float)
                NM_EIGV[mode_IND-1,at,:] = np.array( f[count+1 + at].split()[5:8] ).astype(float)
                NM_EIGV[mode_IND-0,at,:] = np.array( f[count+1 + at].split()[8:11] ).astype(float)

np.save("NM_EIGV.npy", NM_EIGV[:,:,:])

print( good, bad )

if ( Optical_freq is [] ):
    print ("\n\tNo optical frequencies were found. Check calculation.")
    print ("\tSkipping plots.\n")
    exit()

print ( "NModes:", len(FREQ), f" ?= {3*NAtoms} " )
"""
print ("Freq:", FREQ)
print ("IR_INT:", IR_INT)
print ("RAMAN:", RAMAN)
print ("RAMAN Zero Frequency:", RAMAN_f1)
print ("RAMAN First non-Zero Frequency:", RAMAN_f2)
"""
print ( len(FREQ), len(IR_INT), len(RAMAN), len(RAMAN_f1), len(RAMAN_f2) )
print(counter)
np.savetxt("IR_RAMAN_RAW.dat", np.c_[FREQ,IR_INT,RAMAN_f1, RAMAN_f2], header="FREQ, IR, RAMAN (f=0), RAMAN(f!=0)" )

FREQ = np.array(FREQ)
IR_INT = np.array(IR_INT)
RAMAN = np.array(RAMAN)
RAMAN_f2 = np.array(RAMAN_f2)

# Get Spectrum
sig  = 25 # cm^-1, 2.5 meV = 25 cm^-1
Npts = 2000
Emin = 0 # cm^-1
Emax = 2000 # cm^-1
dE = (Emax - Emin) / Npts

IR_SPEC = np.zeros(( Npts ))
RAMAN_SPEC = np.zeros(( Npts ))
RAMAN_f2_SPEC = np.zeros(( Npts ))
DOS = np.zeros(( Npts ))
for j in range( Npts ):
    E = Emin + j*dE
    IR_SPEC[j] += np.sum(  IR_INT * np.exp( -(FREQ - E)**2/2/sig**2 )  )
    RAMAN_SPEC[j] += np.sum(  RAMAN * np.exp( -(FREQ - E)**2/2/sig**2 )  )
    RAMAN_f2_SPEC[j] += np.sum(  RAMAN_f2 * np.exp( -(FREQ - E)**2/2/sig**2 )  )
    DOS[j] += np.sum(  1.000 * np.exp( -(FREQ-E)**2/2/sig**2 )  )


# Save Data to File
np.savetxt( "IR_RAMAN_SPEC.dat", np.array([ np.linspace(Emin,Emax,Npts), IR_SPEC / np.max(IR_SPEC), RAMAN_SPEC / np.max(RAMAN_SPEC), RAMAN_f2_SPEC / np.max(RAMAN_f2_SPEC), DOS / np.max(DOS) ]).T, header="FREQ, IR, RAMAN, RAMAN_f2 (CPHF), DOS" )
np.savetxt( "IR_RAMAN_SPEC_no_NORM.dat", np.array([ np.linspace(Emin,Emax,Npts), IR_SPEC, RAMAN_SPEC, RAMAN_f2_SPEC, DOS ]).T, header="FREQ, IR, RAMAN, RAMAN_f2 (CPHF), DOS" )

plt.plot( np.linspace(Emin,Emax,Npts), IR_SPEC / np.max(IR_SPEC),c="black",label="IR" )
plt.plot( np.linspace(Emin,Emax,Npts), RAMAN_SPEC / np.max(RAMAN_SPEC),c="red",label="RAMAN" )
plt.plot( np.linspace(Emin,Emax,Npts), RAMAN_f2_SPEC / np.max(RAMAN_f2_SPEC),'--',c="red",label=f"CPHF,{round(float(Optical_freq[1])*1.23981*10**-4,2)}eV" )
plt.plot( np.linspace(Emin,Emax,Npts), DOS / np.max(DOS) , linewidth=5, alpha=0.25, c="blue",label="DOS" )
plt.legend()
plt.title(f"Resonant Electronic Transition: {1.2526} eV",fontsize=15)
plt.xlabel("Frequency (1/cm)", fontsize=15)
plt.ylabel("Normalized Intensity", fontsize=15)
plt.savefig("IR_RAMAN_SPEC.jpg")
plt.clf()


# Try plotting in log-scale
plt.plot( np.linspace(Emin,Emax,Npts), np.log(IR_SPEC / np.max(IR_SPEC) + 0.01),c="black",label="IR" )
plt.plot( np.linspace(Emin,Emax,Npts), np.log(RAMAN_SPEC / np.max(RAMAN_SPEC) + 0.01),c="red",label="RAMAN" )
plt.plot( np.linspace(Emin,Emax,Npts), np.log(RAMAN_f2_SPEC / np.max(RAMAN_f2_SPEC) + 0.01),'--',c="red",label="CPHF" )
plt.plot( np.linspace(Emin,Emax,Npts), np.log(DOS / np.max(DOS) + 0.01) , linewidth=5, alpha=0.25, c="blue",label="DOS" )
plt.legend()
plt.xlabel("Frequency (1/cm)", fontsize=15)
plt.ylabel("Log(Intensity)", fontsize=15)
plt.savefig("IR_RAMAN_SPEC_log.jpg")
plt.clf()
















