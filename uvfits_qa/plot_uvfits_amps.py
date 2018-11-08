#!/usr/bin/python

from astropy.io import fits
import sys,os
import numpy as np
from matplotlib import pyplot as plt
import glob

uvlist = glob.glob(sys.argv[1] + '*.uvfits')

# sort by physical frequency

freqs = []
bands = []

for uvfile in uvlist:
    fp = fits.open(uvfile)
    band = uvfile.split('.')[-2][-2:]
    freq = fp[0].header['CRVAL4']
    freqs.append(freq)
    bands.append(band)
    print uvfile, freq, band
    fp.close()

freq_order = np.argsort(freqs)

all_amps = []

for f in freq_order:
    fp = fits.open(uvlist[f])
    data = fp[0].data.data
    #    xx_reals = data[:,0,0,0,:,0,0]
    #    xx_imag = data[:,0,0,0,:,0,1]
    xx_reals = data[:,0,0,:,0,0]
    xx_imag = data[:,0,0,:,0,1]
    xx_weights = data[:,0,0,:,0,2]
    xx = xx_reals.astype(complex) + 1.0j*xx_imag
    xx_amps = np.sqrt(xx_reals * xx_reals + xx_imag * xx_imag)
    xx_amps2 = np.sqrt(xx*np.conjugate(xx))
    flagged_xx_amps = np.ma.masked_where(xx_weights==0,xx_amps)
    #    ch_avg = np.ravel(np.mean(xx_amps,axis=0))
    ch_avg = np.ravel(np.mean(flagged_xx_amps,axis=0))
    all_amps.append(ch_avg)
    fp.close()

all_amps = np.ravel(all_amps)
foo
outfile = 'foo_amps.txt'
out_file = open(outfile,'w+')
for i in all_amps:
    out_file.write('%3.2f ' % i)
out_file.close()
    

plt.clf()
plt.plot(all_amps)
#plt.ylim(max(all_amps) -0.1, max(all_amps) + 0.1)
plt.ylim(10,40)
plt.title(sys.argv[1])
plt.xlabel('Channel Number')
plt.ylabel('Amp(XX)')
#plt.savefig(sys.argv[1] + '.png')
plt.savefig('uvfits_Amps.png')    
    
    
#fp = fits.open(sys.argv[1])

#data = fp[0].data.data

#xx_reals = data[:,0,0,0,:,0,0]
#xx_imag = data[:,0,0,0,:,0,1]

#xx_amps = np.sqrt(xx_reals * xx_reals + xx_imag + xx_imag)

#

#n_baselines = 8128
#n_times = 14

#diff_sum23 = np.zeros(n_times)
#diff_sum34 = np.zeros(n_times)

#for i in range(n_baselines):
#    diff_sum23 += abs(xx_reals[i::8128,2] - xx_reals[i::8128,3]) 
#    diff_sum34 += abs(xx_reals[i::8128,3] - xx_reals[i::8128,4])
    
#xx_reals = np.ravel(xx_reals)
#xx_imag = np.ravel(xx_imag)

#xx_amps = np.sqrt(xx_reals * xx_reals + xx_imag * xx_imag)


