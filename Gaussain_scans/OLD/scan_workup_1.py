# The purpose of this script is to produce a set of 
# single point calcualtions based on scan coordinates
# from Gaussian

# ~ Braden Weight, July 11, 2020  ~ #

# SYNTAX: python3 python_file.py ( geometry.out ) ( TD_Directory_Name )
# (...) Optional argument -- Need to give first if you want to give second but can give first without second

import numpy as np
import re,sys,os,shutil
import subprocess as sp

# Give location of submit script
submit_script="/home/bweight/submit_scripts/submit.gaussian"

# Declare Gaussian output file at command line
if ( len(sys.argv) > 1 ):
  name = sys.argv[1]
else:
  name = "geometry.out"
outFILE = open(name,"r").readlines()

# Make outer directory for single points
if ( len(sys.argv) > 2 ):
  dir1 = f"SP_TDDFT_Singlets_{sys.argv[2]}"
else:
  dir1 = "SP_TDDFT_Singlets_SYMM_x"

QdoDdebug = False

# Read final iterations to relaxed scan steps
file001 = open("Scan_steps_Final.dat","r")
ScanLines = file001.readlines()
iters = []
for i in range(len(ScanLines)):
  iters.append(int(ScanLines[i].split()[2]))
xyz = []
CM = [-1,-1] # Charge and Multiplicity
Natoms = 0
for line in outFILE:
  t = line.split()
  if (len(t) > 0 and t[0] == "NAtoms=" and Natoms == 0):
    Natoms = int(t[1])
    break


counter = 0
cc = 0
symmFlag = False # Checks for symmetry orientation -- 'nosymm' keyword avoids this issue
atomFlag = False # Checks for upcoming geometry
count = 0
for line in outFILE:
  t = line.split()
  if (len(t) > 0 and t[0] == "NAtoms=" and Natoms == 0):
    Natoms = int(t[1])
  if (len(t) > 5 and t[3] == "Multiplicity"):
    CM[0] = t[2]
    CM[1] = t[5]
  if ( len(t) == 2 and t[0] == "Input" and t[1] == "orientation:" ):
    symmFlag = True # Indicates bad orientation
    print (line, symmFlag)
  elif ( len(t) == 2 and t[0] == "Standard" and t[1] == "orientation:" ):
    symmFlag = False # Indicates good orientation
    print (line, symmFlag)
  if (len(t) >= 4 and t[3] == "Coordinates" and symmFlag == True):
    atomFlag = True
    cc = count
  if (len(t) == 6 and atomFlag == True and count - cc >= 3 and count - cc < Natoms + 3):
    at = t[1]
    x = float(t[3])
    y = float(t[4])
    z = float(t[5])
    xyz.append([at,x,y,z])
  if (count - cc == Natoms + 3 and atomFlag == True):
    cc == 0
    atomFlag == False
  count += 1

print ("Natoms", Natoms)
print ("Number of Coordinates Recorded", len(xyz), len(xyz)/Natoms)

# Reshape list to 3D array: [Geoms, Atoms, [Type,x,y,z]]
NGeoms = len(xyz)//Natoms
xyznew = np.array(xyz).reshape(( NGeoms,Natoms,4 ))
sumIters = [ int(np.sum(iters[:it+1]))-1 for it in range(len(iters)) ]
sumIters[-1] = NGeoms-1

#####################
if ( QdoDdebug == True ): # To compare with gaussian, make sure to add "nosymm" to scan .com file
# For debugging, print difference matrix at each step since we are tesing bond scan of At1 and At2
  newScanPoint = []
  for geom in range(1,len(xyznew)):
    print(f"\nGEOM #{geom} - #{geom-1}")
    diffM = np.subtract( np.array(xyznew[geom,:,1:]).astype(float), np.array(xyznew[geom-1,:,1:]).astype(float) )
    print (diffM[:6])
    print (xyznew[geom-1,:2,1:])
    diffR = np.linalg.norm( diffM[29] - diffM[57] )
    print ( "%2.4f" % (diffR) )
    if ( not diffR <= 0.45 ): # Keep this number below scanning interval for this test
      newScanPoint.append(geom-1)

  if ( len(newScanPoint) == 0 or len(newScanPoint) != len(sumIters)-1 ):
    print (f"\tSomething is wrong. Try scanning with 'nosymm' keyword.")
    print ( f"\tLEN(sumIters) {len(sumIters)}, LEN(newScanPoints) {len(newScanPoint)}" )
    exit()

  newScanPoint.append(NGeoms-1)
  print ("New Data", newScanPoint)
  print ("Old Data:",sumIters)

#######################


# Make dictionary for atom types. Modify as needed.
dict = {'6':'C', '1':'H', '7':'N', '8':'O'}

if (not os.path.exists(dir1)):
  os.mkdir(dir1)
else:
  shutil.rmtree(dir1)
  os.mkdir(dir1)

file02 = open('trajectory.xyz','w') # Print all optimized points to file for viewing with VMD or Ovito
for geom in range( len(sumIters) ) :
  print ( f"GEOM: {sumIters[geom]}" )

  dir2 = f"TD_{geom}/"
  totaldir = dir1 + '/' + dir2
  os.mkdir(totaldir)
  file01 = open(totaldir + 'geometry.com','w')
  file01.write("%chk=geometry.chk\n")
  file01.write("%mem=1GB\n")
  file01.write("%nprocshared=1\n\n")

  file01.write("#p CAM-B3LYP/6-31G nosymm\n")
  file01.write("#p TD=(singlets,nstates=15)\n") # Singlet/Triplet Excitations
  file01.write(f"#p IOp(6/8=3) IOp(9/40=4)\n\n") # Print transition density in atomic basis

  file01.write("TitleMe\n\n")

  file01.write("0 1\n") # Charge and multiplicity

  file02.write(str(Natoms) + "\nTitleMe\n")

  for at in range(Natoms):
    line = []
    for k in range(3):
      line.append( "%5.8f" % (np.round( float(xyznew[sumIters[geom],at,k+1]),5) ) ) # Gaussian needs floating, not scientific
    file01.write( dict[str(int(xyznew[sumIters[geom],at,0]))] + "\t\t" + "\t".join(map(str,line)) + '\n' )
    file02.write( dict[str(int(xyznew[sumIters[geom],at,0]))] + "\t\t" + "\t".join(map(str,line)) + '\n' )
  file01.write("\n\n\n\n\n\n\n")
  file01.close()

  # Copy a cluster script to the directory
  sp.call(f"cp {submit_script} {totaldir}", shell=True)
  pwd = os.getcwd()
  os.chdir(totaldir)
  print ("\tSUBMIT")

  sp.call( f"sbatch {submit_script}", shell=True )

  os.chdir(pwd)





















"""



pp1 = -1

file02 = open('trajectory.xyz','w') # Print all optimized points to file for viewing with VMD or Ovito

file03 = open('trajectory_allSteps.xyz','w') # Print all optimized points to file for viewing with VMD or Ovito

for i in range( NGeoms ):

  print (i+1, sum(iters[:pp1+2]) )

  if ((i+1) == sum(iters[:pp1+2]) or i == 0):

    print("\t", newScanPoint[pp1+1], sumIters[pp1+1]) # For debugging

    pp1 += 1

    dir2 = f"TD_{pp1}/"

    totaldir = dir1 + '/' + dir2

    os.mkdir(totaldir)

    file01 = open(totaldir + 'geometry.com','w')

    file01.write("%chk=geometry.chk\n")

    file01.write("%mem=1GB\n")

    file01.write("%nprocshared=1\n\n")

    file01.write("#p B3LYP/STO-3G\n")

    file01.write("#p TD=(singlets,nstates=10)\n") # Singlet/Triplet Excitations

    file01.write(f"#p IOp(6/8=3) IOp(9/40=4)\n\n") # Print transition density in atomic basis

    file01.write("TitleMe\n\n")

    file01.write("0 1\n") # Charge and multiplicity

    #file01.write(CM[0] + ' ' + CM[1] + "\n")

    file02.write(str(Natoms) + "\nTitleMe\n")

    for at in range(Natoms):

      line = []

      for k in range(3):

        line.append( "%5.5f" % (np.round( float(xyznew[i,at,k+1]),5) ) ) # Gaussian needs floating, not scientific

      file01.write( dict[str(int(xyznew[i,at,0]))] + "\t" + "\t".join(map(str,line)) + '\n' )

      file02.write( dict[str(int(xyznew[i,at,0]))] + "\t" + "\t".join(map(str,line)) + '\n' )

    file01.write("\n\n\n\n\n\n\n")

    file01.close()

    # Copy a cluster script to the directory

    sp.call(f"cp {submit_script} {totaldir}", shell=True)

    pwd = os.getcwd()

    os.chdir(totaldir)

    print ("\tSUBMIT")

    #sp.call( f"sbatch {submit_script.split('/')[-1]}", shell=True )

    os.chdir(pwd)



  file03.write(str(Natoms) + "\nTitleMe\n")

  for at in range(Natoms):

    line = []

    for k in range(3):

      line.append( np.round( float(xyznew[i,at,k+1]),5 ) )

    file03.write( dict[str(int(xyznew[i,at,0]))] + "\t" + "\t".join(map(str,line)) + '\n' )



file02.close()

file03.close()



"""
