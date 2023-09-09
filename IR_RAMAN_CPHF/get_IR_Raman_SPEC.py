import numpy as np
from matplotlib import pyplot as plt

# Plotting Parameters
sig  = 25 # cm^-1, 2.5 meV = 25 cm^-1
Npts = 2000
Emin = 0 # cm^-1
Emax = 2000 # cm^-1s


# Read Data
filename = "geometry.out"
f = open(filename,"r") # Open File
FREQ = []
IR_INT = []
RAMAN = []
for line in f:
    t = line.split()
    if ( len(t) == 5 ):
        if ( t[0] == "Frequencies" and t[1] == "--" ):
            FREQ.append( float(t[2]) ) # cm^-1
            FREQ.append( float(t[3]) )
            FREQ.append( float(t[4]) )
    if ( len(t) == 6 ):
        if ( t[0] == "IR" and t[1] == "Inten" ):
            IR_INT.append( float(t[3]) ) # KM/Mole
            IR_INT.append( float(t[4]) )
            IR_INT.append( float(t[5]) )
        if ( t[0] == "Raman" and t[1] == "Activ" ):
            RAMAN.append( float(t[3]) ) # A**4 / AMU
            RAMAN.append( float(t[4]) )
            RAMAN.append( float(t[5]) )        

print ( "NModes:", len(FREQ), " = 3*NAtoms ?" )
print ("Freq:", FREQ)
print ("IR_INT:", IR_INT)
print ("RAMAN:", RAMAN)

FREQ = np.array(FREQ)
IR_INT = np.array(IR_INT)
RAMAN = np.array(RAMAN)

# Get Spectrum
dE = (Emax - Emin) / Npts

IR_SPEC = np.zeros(( Npts ))
RAMAN_SPEC = np.zeros(( Npts ))
DOS = np.zeros(( Npts ))
for j in range( Npts ):
    E = Emin + j*dE
    IR_SPEC[j] += np.sum(  IR_INT * np.exp( -(FREQ - E)**2/2/sig**2 )  )
    if (len(RAMAN) > 0): RAMAN_SPEC[j] += np.sum(  RAMAN * np.exp( -(FREQ - E)**2/2/sig**2 )  )
    DOS[j] += np.sum(  1.000 * np.exp( -(FREQ-E)**2/2/sig**2 )  )

plt.plot( np.linspace(Emin,Emax,Npts), IR_SPEC / np.max(IR_SPEC),c="black",label="IR" )
if (len(RAMAN) > 0): plt.plot( np.linspace(Emin,Emax,Npts), RAMAN_SPEC / np.max(RAMAN_SPEC),c="red",label="RAMAN" )
plt.plot( np.linspace(Emin,Emax,Npts), DOS / np.max(DOS) , "--",c="blue",label="DOS" )
plt.legend()
plt.xlabel("Frequency (1/cm)", fontsize=15)
plt.ylabel("Normalized Intensity", fontsize=15)
plt.savefig("IR_RAMAN_SPEC.jpg")
plt.clf()

# Try plotting in log-scale
plt.plot( np.linspace(Emin,Emax,Npts), np.log(IR_SPEC / np.max(IR_SPEC) + 0.01),c="black",label="IR" )
if (len(RAMAN) > 0): plt.plot( np.linspace(Emin,Emax,Npts), np.log(RAMAN_SPEC / np.max(RAMAN_SPEC) + 0.01),c="red",label="RAMAN" )
plt.plot( np.linspace(Emin,Emax,Npts), np.log(DOS / np.max(DOS) + 0.01) , "--",c="blue",label="DOS" )
plt.legend()
plt.xlabel("Frequency (1/cm)", fontsize=15)
plt.ylabel("Log(Intensity)", fontsize=15)
plt.savefig("IR_RAMAN_SPEC_log.jpg")
plt.clf()

if (len(RAMAN) > 0):
    np.savetxt("FREQ_IR_RAMAN.dat", np.c_[FREQ, IR_INT, RAMAN] )
    np.savetxt("FREQ_DOS_IR_RAMAN_SPEC.dat", np.c_[np.linspace(Emin,Emax,Npts), DOS, IR_SPEC, RAMAN_SPEC] )
else:
    np.savetxt("FREQ_IR.dat", np.c_[FREQ, IR_INT] )
    np.savetxt("FREQ_DOS_IR_SPEC.dat", np.c_[np.linspace(Emin,Emax,Npts), DOS, IR_SPEC] )














