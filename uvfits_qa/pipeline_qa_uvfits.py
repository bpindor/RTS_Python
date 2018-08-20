#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 14:06:19 2014

@author: bpindor
Reads values from RTS uvfits files
"""

from astropy.io import fits
import sys,os
from numpy import zeros, shape, real, imag, arange, ravel, max, sqrt, std, empty, concatenate, where, argsort

def pipeline_qa_uvfits(obsid,options):

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

    n_bands = options.nbands
    n_channels = 16 # 32 / 2, flagged channels present
    sigmas = empty((n_bands,4))
    ch_sigmas = []

    freq_list = []

    # sort bands by frequency

    for band in range(n_bands):

        infile = uvfits_dir + 'uvdump_%02d.uvfits' % (band+1)

        fp = fits.open(infile)

        freq_list.append(fp[0].header['CRVAL4'])
    
        fp.close()

    band_order = argsort(freq_list)   
    
    for band in range(n_bands):

        if(options.tagname is None):
            infile = uvfits_dir + ('uvdump_%02d' % (band_order[band]+1)) + '.uvfits'
        else:
            infile = uvfits_dir + ('%s_%02d' % (options.tagname,band_order[band]+1)) + '.uvfits'

        print infile

        fp = fits.open(infile)

        n_groups = fp[0].header['GCOUNT']
        n_ant = fp[1].header['NAXIS2']
        n_chan = fp[0].header['NAXIS4']
        n_baselines = n_ant * (n_ant - 1) / 2
# This should be an option or better still somehow calculated from metadata
        n_time = 14 # Number of (8s) timesteps
        n_diffs = n_time / 2 # How many differences can we form

        xx_reals = []
        xx_imag = []
        yy_reals = []
        yy_imag = []
        weights = []

        data = fp[0].data.data

        xx_reals = data[:,0,0,:,0,0]
        xx_imag = data[:,0,0,:,0,1]
        yy_reals = data[:,0,0,:,3,0]
        yy_imag = data[:,0,0,:,3,1]
        weights = data[:,0,0,:,0,2]

        xx_reals = ravel(xx_reals)
        xx_imag = ravel(xx_imag)
        yy_reals = ravel(yy_reals)
        yy_imag = ravel(yy_imag)
        weights = ravel(weights)
                
        # Maybe this differencing could be done by reading two sets of data?

        for i in range(n_diffs):
            if(i==0):
                xxr_diffs = xx_reals[i*2*n_baselines*n_chan:(i*2+1)*n_baselines*n_chan] - xx_reals[(i*2+1)*n_baselines*n_chan:(i*2+2)*n_baselines*n_chan]
                xxi_diffs = xx_imag[i*2*n_baselines*n_chan:(i*2+1)*n_baselines*n_chan] - xx_imag[(i*2+1)*n_baselines*n_chan:(i*2+2)*n_baselines*n_chan]
                yyr_diffs = yy_reals[i*2*n_baselines*n_chan:(i*2+1)*n_baselines*n_chan] - yy_reals[(i*2+1)*n_baselines*n_chan:(i*2+2)*n_baselines*n_chan]
                yyi_diffs = yy_imag[i*2*n_baselines*n_chan:(i*2+1)*n_baselines*n_chan] - yy_imag[(i*2+1)*n_baselines*n_chan:(i*2+2)*n_baselines*n_chan] 
                weights_list = sqrt(weights[i*2*n_baselines*n_chan:(i*2+1)*n_baselines*n_chan] * weights[(i*2+1)*n_baselines*n_chan:(i*2+2)*n_baselines*n_chan])
            else:
                xxr_diffs = concatenate((xxr_diffs,xx_reals[i*2*n_baselines*n_chan:(i*2+1)*n_baselines*n_chan] - xx_reals[(i*2+1)*n_baselines*n_chan:(i*2+2)*n_baselines*n_chan]))
                xxi_diffs = concatenate((xxi_diffs,xx_imag[i*2*n_baselines*n_chan:(i*2+1)*n_baselines*n_chan] - xx_imag[(i*2+1)*n_baselines*n_chan:(i*2+2)*n_baselines*n_chan]))                        
                yyr_diffs = concatenate((yyr_diffs,yy_reals[i*2*n_baselines*n_chan:(i*2+1)*n_baselines*n_chan] - yy_reals[(i*2+1)*n_baselines*n_chan:(i*2+2)*n_baselines*n_chan]))
                yyi_diffs = concatenate((yyi_diffs,yy_imag[i*2*n_baselines*n_chan:(i*2+1)*n_baselines*n_chan] - yy_imag[(i*2+1)*n_baselines*n_chan:(i*2+2)*n_baselines*n_chan]))                        
                weights_list = concatenate((weights_list,sqrt(weights[i*2*n_baselines*n_chan:(i*2+1)*n_baselines*n_chan] * weights[(i*2+1)*n_baselines*n_chan:(i*2+2)*n_baselines*n_chan]))) 

        # Array now contain visibility differences, just need to calculate variances

        rxx = xxr_diffs
        ixx = xxi_diffs
        ryy = yyr_diffs
        iyy = yyi_diffs
        
        
        r_weights = ravel(weights_list)
        max_weight = r_weights.max()
        nz_weights = where(r_weights > 0)

        weighted_xxr = xxr_diffs[nz_weights] / sqrt(max_weight / r_weights[nz_weights])
        weighted_xxi = xxi_diffs[nz_weights] / sqrt(max_weight / r_weights[nz_weights])
        weighted_yyr = yyr_diffs[nz_weights] / sqrt(max_weight / r_weights[nz_weights])
        weighted_yyi = yyi_diffs[nz_weights] / sqrt(max_weight / r_weights[nz_weights])

        sigmas[band] = [std(weighted_xxr),std(weighted_xxi),std(weighted_yyr),std(weighted_yyi)]
        for n in range(n_channels):
            xxr_ch = rxx[n::16]
            weights_ch = r_weights[n::16]
            nz_weights = where(weights_ch > 0)
            if(len(nz_weights[0]) < 1):
                ch_sigmas.append(0.0)
            else:
                #weighted_xxr = xxr_ch[nz_weights]
                weighted_xxr = xxr_ch[nz_weights] / sqrt(max_weight / weights_ch[nz_weights])
                ch_sigmas.append(std(weighted_xxr))
                
#        sigmas.append([std(weighted_xxr),std(weighted_xxi),std(weighted_yyr),std(weighted_yyi)])

        fp.close()

    outfile = obs_num + '_rms.txt'

    out_file = open(outfile,'w+')

    ch_sigmas = ravel(ch_sigmas)


    for i in range(len(ch_sigmas)):
       out_file.write('%3.2f ' % ch_sigmas[i]) 
    
#    for i in range(n_coarse):
#        out_file.write('%3.2f ' % sigmas[i,0])


    out_file.close()



        
    



from optparse import OptionParser,OptionGroup

mwa_dir = os.getenv('MWA_DIR','/scratch2/mwaeor/MWA/')

usage = 'Usage: pipeline_qa.py <obsid_list>'

parser = OptionParser(usage=usage)

parser.add_option('--tagname', dest='tagname',default=None,
                      help='Tag string used to identify the uvfits file to be processed [default=%default]')

parser.add_option('--subdir',dest="subdir",type='string',default='',
                  help="Data subdirector where RTS will be run and outputs stored")
parser.add_option('--nbands',dest="nbands",type='int',default='24',
                  help="Number of Coarse Channels Present")

(options, args) = parser.parse_args()

obsid = args[0]

pipeline_qa_uvfits(obsid,options)
