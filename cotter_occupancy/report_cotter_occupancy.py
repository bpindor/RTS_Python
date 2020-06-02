from __future__ import print_function
from astropy.io import fits 
import sys,os
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from optparse import OptionParser,OptionGroup
import numpy as np

usage = 'Usage: scan_mwaf_files.py [text file of obsIDs]'

parser = OptionParser(usage=usage)

parser.add_option('-t','--threshold',dest='threshold',default=0.8,type='float',help='Fraction of flagged channels which triggers flagging of entire channel')
parser.add_option('-d','--dir',dest='mwa_dir',default='',type='string',help='Base Data Directory')
parser.add_option('-n','--n_bands',dest='n_bands',default=24,type='int',help='Number of Coarse Bands to Use')
parser.add_option('-p','--prefix',dest='prefix',default=None,type='string',help='MWAF files prefix string')
parser.add_option('-v','--verbose',dest='verbose',action='store_true',default=False,help='MWAF files prefix string')


(options, args) = parser.parse_args()

obsdir = args[0]

obsid = obsdir.split('/')[0]

if options.mwa_dir:
    mwa_dir = options.mwa_dir
else:
    mwa_dir = os.getenv('MWA_DIR','/scratch2/astronomy883/MWA/')

data_dir = mwa_dir + 'data/%s/' % obsid

print(data_dir)

if(options.prefix is None):
    prefix = '1'
    outfile = 'CotterOccupancy_%s.dat' % obsid
    plotfile = 'CotterOccupancy_%s.png' % obsid
else:
    prefix = options.prefix + '_1'
    outfile = options.prefix + '_CotterOccupancy_%s.dat' % obsid
    plotfile = options.prefix + '_CotterOccupancy_%s.png' % obsid

mwaf_list = glob.glob(data_dir + '/%s*.mwaf' % prefix)

# remove flagged channels
channels = np.arange(2,30)
channels = np.delete(channels,14)

if(len(mwaf_list) != 24):
    print('Error: %s had %d mwaf files' % (obsid, len(mwaf_list)))
    exit(1)

out_file = open(outfile,'w+')    
    
unflagged_baselines = None
unflagged_timesteps = None

mwaf_list.sort()

#plt.clf()
plt.subplot2grid((2,1),(0,0))

all_occupancy = []

for mwaf_file in mwaf_list: 

    if(options.verbose):
        print(mwaf_file)
    band = (mwaf_file.split('_')[-1])[0:2]
    if(int(band) > options.n_bands):
        continue
    out_file.write('%s ' % band)

    hdulist = fits.open(mwaf_file)

    flags=hdulist[1].data['FLAGS']
    n_ants = hdulist[0].header['NANTENNA']

    # sometimes the actual number of scans is less than n_scans
    n_scans = hdulist[0].header['NSCANS']
    n_baselines = (n_ants * (n_ants + 1)) / 2
    n_scans = int((flags.shape[0]) / n_baselines)

    if(options.verbose):
        print('NSCANS',n_scans)

    n_chan = np.shape(flags)[1]

    channels = np.arange(n_chan/16, n_chan - n_chan/16,dtype=int)
    channels = np.delete(channels,int(n_chan/2 - n_chan/16))

    if(unflagged_baselines == None):
        unflagged_baselines = [i for i in range(8256) if (np.sum(flags[i::8256,channels])) < 0.9 * float(len(channels) * n_scans)]

    flags_per_timestep = [np.sum(flags[n*8256:(n+1)*8256,channels]) for n in range(n_scans)]
    unflagged_timesteps = [n for n in range(n_scans) if flags_per_timestep[n] < 8256 * len(channels) * 0.9]

    unflagged_indices = [np.array(unflagged_baselines) + n * 8256 for n in unflagged_timesteps]

    unflagged_indices = np.ravel(unflagged_indices)

    if(options.verbose):
        print('Unflagged indices: ',len(unflagged_indices))

    # Cant seem to slice on baseline index and channels at the same time
    # So do one after the other
    
    cut_1 = flags[unflagged_indices,:]
    cut_flags = cut_1[:,channels]
    
#    cut_flags = flags[unflagged_timesteps[0]*8256:(unflagged_timesteps[-1]+1)*8256,channels]

    n_flags = []

    n_ch_flags = []

    threshold = float(options.threshold)

    #n_ch_flags = np.sum(flags,axis=0)
    n_ch_flags = np.sum(cut_flags,axis=0)

#    for i in range(32):
#        ch_flags = [f[i] for f in flags]
#        n_ch_flags.append(sum(ch_flags))

    for i in range(len(n_ch_flags)):
        out_file.write('%d ' % n_ch_flags[i])
    
    out_file.write('\n')
    
#    occupancy1 = [float(n) / float(n_ch_flags[0]) for n in n_ch_flags[channels] ]
    occupancy1 = [float(n) for n in n_ch_flags]
    
#    plt.ylim(0,1.1)
    #plt.bar(range(27),occupancy1)
    plt.plot(occupancy1)
    all_occupancy.append(occupancy1)

plt.title('Cotter Occupancy %s ' % obsid)    
plt.subplot2grid((2,1),(1,0))
plt.ylabel('Percent Flagged')
plt.plot(100.0 * (np.ravel(all_occupancy)) / float((len(unflagged_indices))))
plt.savefig('%s' % plotfile)
out_file.close()

percent_file = open('Cotter_percentage_%s.dat' % obsid,'w+')
print(100.0 * (np.ravel(all_occupancy)) / float((len(unflagged_indices))),file=percent_file)
percent_file.close()

flags_per_baseline = [np.sum(flags[i::8256,channels]) for i in unflagged_baselines]
flags_per_timestep = [np.sum(flags[n*8256:(n+1)*8256,channels]) for n in range(n_scans)]

