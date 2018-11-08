#!/usr/bin/python

from astropy.io import fits
import sys,os
import numpy as np
from matplotlib import pyplot as plt
import glob

def pipeline_uvfits_amps(obsid,options):

    #define directories
    mwa_dir = os.getenv('MWA_DIR','/scratch/partner1019/MWA/')
    working_dir= os.getcwd() + "/"
    obs_num=str(obsid)

    if(options.tagname is None):
        if(options.subdir is None):
            uvfits_dir = mwa_dir + 'data/' + obs_num + '/'
        else:
            uvfits_dir = mwa_dir + 'data/' + obs_num + '/%s/' % options.subdir 
    else:
        uvfits_dir = mwa_dir + 'data/' + obs_num + '/uvfits/'
    
    data_dir=mwa_dir+ 'data/'
    

    uvlist = glob.glob(uvfits_dir + '*.uvfits')

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
        max_weight = xx_weights.max()
        xx = xx_reals.astype(complex) + 1.0j*xx_imag
        xx_amps = np.sqrt(xx_reals * xx_reals + xx_imag * xx_imag)
        xx_amps2 = np.sqrt(xx*np.conjugate(xx))
        flagged_xx_amps = np.ma.masked_where(xx_weights==0,xx_amps)
        weighted_xx_amps = flagged_xx_amps / np.sqrt(xx_weights / max_weights)
        #    ch_avg = np.ravel(np.mean(xx_amps,axis=0))
        #ch_avg = np.ravel(np.mean(flagged_xx_amps,axis=0))
        ch_avg = np.ravel(np.mean(weighted_xx_amps,axis=0)) 
        all_amps.append(ch_avg)
        fp.close()

    all_amps = np.ravel(all_amps)

    outfile = '%s_amps.txt' % obsid
    out_file = open(outfile,'w+')
    for i in all_amps:
        out_file.write('%3.2f ' % i)
    out_file.close()
    

#plt.clf()
#plt.plot(all_amps)
#plt.ylim(max(all_amps) -0.1, max(all_amps) + 0.1)
#plt.ylim(10,40)
#plt.title(sys.argv[1])
#plt.xlabel('Channel Number')
#plt.ylabel('Amp(XX)')
#plt.savefig(sys.argv[1] + '.png')
#plt.savefig('uvfits_Amps.png')    

##########################################

from optparse import OptionParser,OptionGroup

mwa_dir = os.getenv('MWA_DIR','/scratch2/mwaeor/MWA/')

usage = 'Usage: pipeline_uvfits_amps.py <obsid_list>'

parser = OptionParser(usage=usage)

parser.add_option('--tagname', dest='tagname',default=None,
                      help='Tag string used to identify the uvfits file to be processed [default=%default]')

parser.add_option('--subdir',dest="subdir",type='string',default='',
                  help="Data subdirector where RTS will be run and outputs stored")
parser.add_option('--nbands',dest="nbands",type='int',default='24',
                  help="Number of Coarse Channels Present")

(options, args) = parser.parse_args()

obsid = args[0]

pipeline_uvfits_amps(obsid,options)

