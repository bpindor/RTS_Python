import RTS_cals
from astropy.io import fits
import os,sys
from optparse import OptionParser,OptionGroup
import numpy as np
from matplotlib import pyplot as plt

usage = 'Usage: single_rewrite_cals [directory containing RTS cals] [metafits]'

parser = OptionParser(usage=usage)

parser.add_option('--n_bands',dest='n_bands',type='int',default=24,help='Number of Coarse Bands Present')
parser.add_option('--plot_tile',dest='plot_tile',type='int',default=-1,help='Number of RTS Tile to Plot')
parser.add_option('--limits',dest='plot_limits',type=string

(options, args) = parser.parse_args()

cal_dir = args[0]
metafile = args[1]

n_bands = options.n_bands

rts_cal = RTS_cals.rts_cal(n_bands=n_bands)
rts_cal.load_all_BP_jones(path=cal_dir, raw=True,n_bands=n_bands)
rts_cal.load_all_DI_jones(path=cal_dir,n_bands=n_bands)
# Form product of RTS BP and DI cals to get a single Jones matrix per tile
# per fine channel
rts_cal.form_single_jones()

meta_file = fits.open(metafile)
cent_chan = meta_file[0].header['CENTCHAN']
flags = meta_file[1].data['Flag']
flags = flags[::2]
unflagged = [i for i in range(rts_cal.n_ants) if i not in flags]

all_single_jones = [None] * rts_cal.n_ants

# The cals are still per coarse in coarse band ordering. Now we form a single
# array for each tile
for i in range(rts_cal.n_ants):
    if i not in flags:
        single_jones = rts_cal.all_single_jones(antenna=i,cent_chan=cent_chan)
        all_single_jones[i] = single_jones
        
# Now each unflagged tile has a full band calibration with dimensions (2,2,# of fine channels).
# Plot real and imaginary of XX for first unflagged tile (by default)

if(options.plot_tile == -1):
    plot_index = unflagged[0]
else:
    plot_index = options.plot_tile
    if(plot_index in flagged):
        print 'ERROR: Cannot plot solutions of flagged tile'
        exit()

plt.clf()
plt.plot(np.real(all_single_jones[plot_index][0][0][:]))
plt.plot(np.imag(all_single_jones[plot_index][0][0][:]))
plt.title('RTS Tile %d' % plot_index)
plt.ylim(0.95,1.05)
plt.savefig('rts_cal.png')
        
