import numpy as np
import imageio

METHOD = 'better' # 'fast' or 'better' # Use 'better' for long videos
filenames = [ f"./NAME_{j}.jpg" for j in range() ] # Change this as needed
movie_Name = './' + 'myNewGIF.gif' # Change this as needed



if ( METHOD == 'fast' ):
    # Fast, easy version  
    imageio.mimsave(movie_Name, filenames)
elif ( METHOD == 'better' ):
    with imageio.get_writer(movie_Name, mode='I', fps=0.1) as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)