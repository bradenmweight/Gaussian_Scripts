# The purpose of this script is to read-in a set of JPEG files
#   and export a video in GIF format.

# Looking in the paths: TD_*/TDMat/trans-*.JPEG

####### Braden M. Weight ~ July 12, 2021 ######

"""
Install pygifcicle:
pip3 install pygifsicle


Install gifsicle: ( echo "$(pwd)" = /scratch/bweight/software/ )
curl -sL http://www.lcdf.org/gifsicle/gifsicle-1.91.tar.gz | tar -zx
cd gifsicle-1.91
./configure --disable-gifview
make install exec_prefix=$(pwd) prefix=$(pwd) datarootdir=$(pwd)
"""

import numpy as np
import imageio
import os
import subprocess as sp
from pygifsicle import optimize as gifOPT # This needs to be installed somewhere

NSteps = 750 # Number of TD_* directories
movieName = "GIF" # "GIF" --> ./GIF.gif
PATH_TO_VMD_CUBE_GEN="~/Gaussain_scripts/Cube_Images/run_cubes_vmd_cubes.sh" # This file makes JPEGS for all *cube files in current directory -- Very Expensive, best to do beforehand

TDDirs = [ f"TD_{j}" for j in range(NSteps) ] # Compile outer directory list
activeStates = [ open(f"TD_{j}/Active_State.dat","r").readlines()[1].split()[1] for j in range(NSteps) ] # Get active state in each directory
filenames = [ f"{TDDirs[j]}/TDMat/trans-{activeStates[j]}.jpg" for j in range(NSteps) ] # Get file paths for the appropriate transition density JPEG file
moviePATH = f'./{movieName}.gif' # Set relative path to move file


with imageio.get_writer(moviePATH, mode='I', fps=5) as writer: # Get a writer object
    for filename in filenames: # Loop over JPEGs
        print ("Filename:", filename)
        if ( not os.path.isfile(filename) ): # Check whether we want to make new JPG files or not. If it doesn't exist, it will error out.
            cwd = os.path.cwd()
            os.chdir( "/".join(filename.split("/")[:-1]) )
            sp.call(PATH_TO_VMD_CUBE_GEN,shell=True) # Change this file to edit the isovalue. For resolution, change ~/.../vmd_cube.py
            os.chdir(cwd)
        image = imageio.imread(filename) # Read JPEG file
        writer.append_data(image) # Write JPEG file (to memory at first; then printed at end)

gifOPT(moviePATH) # This will compress the GIF movie by at least a factor of two/three.