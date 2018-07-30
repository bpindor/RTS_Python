import sys
from mwa_baseline_conversions import *
from astropy.io import fits
import numpy as np

# Reads in a single uvfits file, flags one of the fine (80khz) channels, and writes out a new copy of the data

uvfits = fits.open(sys.argv[1])

# 
uvdata = uvfits[0].data.data

print np.shape(uvdata)
print '(Number of Baselines (8128) * Number of TimeSteps (14), 1, 1, Number of Frequency Channels, Number of Polarizations, [Real, Imag, Weight])'

weights = uvdata[:,0,0,:,0,2]

print np.shape(weights)

# Tricky bit that uvfits are stored by baseline, so we want to flag all of the baselines corresponding to a given tile

# This functions load a bunch of baseline ordering functions based on the metafits

metafits = sys.argv[2]
loadCotterMapping(metafits)

# Suppose we want to flag tile104

print tile2rts(104)

# Tile104 is 0 in rts order (which is order that uvfits are written)

baselines = rtsbaselinesForAnt(tile2rts(104))
baselines = np.array(baselines)

# Now manually set those weights to zero

# Slightly-kludgy loop over timesteps

n_times = 14
n_baselines = 8128

all_baselines = [b + n * n_baselines for b in baselines for n in range(n_times)]

# Suppose we want to flag channel 4

weights[all_baselines,4] = 0.0

uvdata[:,0,0,:,0,2] = weights

# Now we just need to write the data back out

#uvfits[0].data.data = uvdata

uvfits.writeto('flagged.uvfits',clobber=True)

uvfits.close()


    

