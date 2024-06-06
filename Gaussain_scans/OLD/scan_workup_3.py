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
import os, sys
import subprocess as sp
import multiprocessing as mp
from PIL import Image, ImageDraw, ImageFont
from pygifsicle import optimize as gifOPT # This needs to be installed somewhere

NSteps = 20 # Number of TD_* directories
NStates = 10 # This will look for the lowest {NStates} excited states and make three movies, one for each
NSkip = 1 # Compile every {NSkip} frames
FPS = 2 # Frames per second ~ 2-15
movieName = "SCAN_GIF" # "GIF" --> ./GIF_S{n}.gif, n is label of state
PATH_TO_VMD_CUBE_GEN="~/Gaussain_scripts/Cube_Images/run_cubes_vmd_cubes.sh" # This file makes JPEGS for all *cube files in current directory -- Very Expensive
FontFile="/scratch/bweight/anaconda3/lib/python3.8/site-packages/anaconda_navigator/static/fonts/UbuntuMono-R.ttf"


TDDirs = [ f"TD_{j}" for j in range(NSteps) ] # Compile outer directory list
filenames = [ f"{TDDirs[j]}/TDMat/trans-{state+1}.jpg" for j in range(NSteps) for state in range(NStates) ] # Get file paths for the appropriate transition density JPEG file


def addText(fileloc,text):
    img = Image.open(fileloc)
    d1 = ImageDraw.Draw(img)
    d1.text( (10, 0), text, (0,0,0), font=ImageFont.truetype(FontFile,50) )
    img.save(fileloc)


def makeJPGs(filename): # This function is parallelized (optional)
    if ( not os.path.isfile(filename) ): # Check whether we want to make new JPG files or not.
        print (f"\tMaking JPEG in {filename}")
        cwd = os.getcwd()
        os.chdir( "/".join(filename.split("/")[:-1]) )
        sp.call(PATH_TO_VMD_CUBE_GEN,shell=True) # Change this file to edit the isovalue. For resolution, change ~/.../vmd_cube.py
        addText( filename.split("/")[-1], f" Step: {filename.split('/')[0].split('_')[1]}" )
        os.chdir(cwd)
    else:
        print (f"\t\t{filename} is done.")

# Make JPEGs -- This is expensive -- Option for multi-threading by setting NCPUS >= 2
if ( len(sys.argv) == 2 ):
    NCPUS = int(sys.argv[1])
else:
    NCPUS = 1 # Default is single core

with mp.Pool(processes=NCPUS) as pool:
    pool.map(makeJPGs,filenames)


for state in range(NStates):

    moviePATH = f'./{movieName}_S{state}.gif' # Set relative path to move file
    
    # Compile all JPEGS -- This is cheaper
    with imageio.get_writer(moviePATH, mode='I', fps=FPS) as writer: # Get a writer object
        for filename in filenames[state::NStates]: # Loop over JPEGs, NSkip will only keep every 'NSkip' frames
            print ("Filename:", filename)
            try:
                image = imageio.imread(filename) # Read JPEG file
            except ValueError:
                print ("!!!!! Bad JPEG file. Trying to remake it. !!!!!")
                cwd = os.getcwd()
                os.chdir( "/".join(filename.split("/")[:-1]) )
                sp.call(PATH_TO_VMD_CUBE_GEN,shell=True) # Change this file to edit the isovalue. For resolution, change ~/.../vmd_cube.py
                addText( filename.split("/")[-1], f" Step: {filename.split('/')[0].split('_')[1]}" )
                os.chdir(cwd)
                image = imageio.imread(filename) # Read JPEG file
            
            writer.append_data(image) # Write JPEG file (to memory at first; then printed at end)

    gifOPT(moviePATH) # This will compress the GIF movie by at least a factor of two/three. With this: ~750 frames --> 80 MB