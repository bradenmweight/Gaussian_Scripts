import numpy as np

data     = np.loadtxt("IR_RAMAN_RAW.dat")

data = data[ data[:, 3].argsort() ]

np.savetxt("IR_RAMAN_RAW_f1SORTED.dat", data)