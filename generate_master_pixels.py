import numpy as np
from pathlib import Path
import read_params

def populate(num_src):
	datadir = Path(read_params.get_directory())
	# Choose uniformly random numbers between [-20,20)
	master_pixels = np.random.random((num_src,1))*40 - 20
	
	print("master pixels",master_pixels.flatten())
	print("Writing to {}".format(datadir/"master.pixels"))
	
	np.savetxt(datadir/"master.pixels",master_pixels,fmt="%.1f")

