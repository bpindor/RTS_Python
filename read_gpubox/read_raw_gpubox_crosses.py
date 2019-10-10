from matplotlib import pyplot as plt
from load_gpubox_vis import *
import numpy as np
from glob import glob
from spoof_mwac import *
from astropy.io import fits
from fit_cable_reflections import *
from mwa_baseline_conversions import *
from load_mwaf_flags import *

n_iter= 32
nchan = 32
n_ants = 128
n_bands = 20

# This creates a mapping between the rts baselines (which are indexed according to the metafits 'INPUT' value)

corr_mapping = fill_mapping_matrix()
the_matrix = extract_full_matrix()
rts2index = load_rts_baselines(corr_mapping,the_matrix)
#rts_baselines = range(n_ants * (n_ants - 1) / 2)

#rts_baselines = range(n_ants * (n_ants - 1) / 2)

rts_baselines = range(n_ants * (n_ants - 1) / 2)
raw_baselines = [rts2index[i] for i in rts_baselines]

# Creates a conversion between the Cotter baselines (used in the mwaf files) and the rts baselines

#loadCotterMapping(meta_file)

#cotter_baselines = [rtsbl2Cotter(i) for i in rts_baselines]

#

mwa_dir = os.getenv('MWA_DIR','/astro/mwaeor/MWA/')

obsid = sys.argv[1]

gpu_files = glob(mwa_dir+'data/'+obsid+'/*gpubox*00.fits')
mwaf_files = glob(mwa_dir+'data/'+obsid+'/%s*.mwaf' % obsid)

metafiles = glob(mwa_dir+'data/'+obsid+'/*meta*fits')
meta_file = fits.open(metafiles[0])
cent_chan = meta_file[0].header['CENTCHAN']
band_0 = cent_chan - 12

loadCotterMapping(metafiles[0])

cotter_baselines = [rtsbl2Cotter(i) for i in rts_baselines]


rts_inputs = meta_file[1].data['Input']
metafiles = glob(mwa_dir+'data/'+obsid+'/*meta*ppds*fits')
meta_file = fits.open(metafiles[0])
cables = meta_file[1].data['Flavors']
cable_types = set(cables)
#cable_counts = dict(Counter(cables[::2])) 
rts_ants = rts_inputs[::2] / 2
rts2cables = dict(zip(rts_ants,cables[::2]))
flags = meta_file[1].data['Flag']
flags = flags[::2]
channels = meta_file[0].header['Channels']
channels = channels.split(',')


# Digital gains

cable2gains = {}

gain_flavors = meta_file[2].data['Flavors']
digital_gains = meta_file[2].data['Gains'] 

for i,flavor in enumerate(gain_flavors):
    cable2gains[flavor] = digital_gains[i][int(channels[0]):int(channels[0])+n_bands]

gpu_files.sort()
mwaf_files.sort()

band_numbers = [int((g.split('gpubox')[1])[:2])-1 for g in gpu_files]

gpu_files = gpu_files[:n_bands]

# Here we start reading the data

#flag_vis = zeros((n_iter,nchan*n_bands,len(rts_baselines)),dtype=complex)
raw_vis = zeros((n_iter,nchan*n_bands,len(rts_baselines)),dtype=complex)

#for band,gpubox_file in enumerate(gpu_files):
for band,gpubox_file in zip(band_numbers,gpu_files):
    
    all_vis = get_raw_vis(gpubox_file,n_iter=n_iter)

    all_flags = load_mwaf_flags(mwaf_files[band])

    print gpubox_file,mwaf_files[band]

    # Correct band ordering
    
    band_number = band + cent_chan - 12
    if(band_number > 128):
        if(band_0 > 128):
            band_index = n_bands -(band+1)
        else:
            band_index = n_bands-(band_number-128)
    else:
        band_index = band
     
    print band_number, band_index

    dg_term = np.zeros(len(rts_baselines))
    
    for b in rts_baselines:
        (ant1,ant2) = rtsbl2rtsantennas[b]
        dg1 = cable2gains[rts2cables[ant1]][band_index]
        dg2 = cable2gains[rts2cables[ant2]][band_index]
        dg_term[b] = ((float) (dg1 * dg2) / (64.0)**2)
    
    for t in range(n_iter):
        for ch in range(nchan):
            for b in rts_baselines:
                raw_vis[t,band_index*nchan+ch,b] = all_vis[t][ch][raw_baselines[b]]
#                raw_vis[t,band_index*nchan+ch,b] /= dg_term[b]
                if(all_flags[t,cotter_baselines[b],ch] == True):
                    raw_vis[t,band_index*nchan+ch,b] = 0.0
                
mean_vis = np.mean(abs(raw_vis),axis=0)

#bl_flags = np.zeros(len(rts_baselines))

low_band_pfb = np.array([ 0.52200223,  0.58254809,  0.68921804,  0.80630983,  0.90282802,
        0.96484245,  0.99443162,  1.00278847,  1.00158203,  0.9993008 ,
        0.99899039,  1.00050617,  1.00177455,  1.00171228,  1.00036634,
        0.99899723,  1.00135908,  0.9997995 ,  1.00112694,  1.00141863,
        1.00038367,  0.99874683,  0.99849967,  1.00019452,  1.00232282,
        0.99707395,  0.97391232,  0.9206596 ,  0.83210091,  0.71841159,
        0.60662328,  0.53326278])

#pfb_corrected_vis = apply_pfb_to_autos(raw_vis,bl_flags,low_band_pfb,n_bands=n_bands)

#mean_corrected_vis = np.mean(abs(pfb_corrected_vis),axis=0)

emp_pfb = np.array([0.5       ,  0.5       ,  0.67874855,  0.83576969,  0.95187049,
        1.0229769 ,  1.05711736,  1.06407012,  1.06311151,  1.06089592,
        1.0593481 ,  1.06025714,  1.06110822,  1.05893943,  1.05765503,
        1.05601938,  0.5       ,  1.05697461,  1.05691842,  1.05688129,
        1.05623901,  1.05272663,  1.05272112,  1.05551337,  1.05724941,
        1.0519857 ,  1.02483081,  0.96454596,  0.86071928,  0.71382954,
					    0.5       ,  0.5])

#emp_corrected_vis = apply_pfb_to_autos(raw_vis,bl_flags,emp_pfb,n_bands=n_bands)

#mean_emp_corrected_vis = np.mean(abs(emp_corrected_vis),axis=0)

#cross_diffs =  np.zeros((n_iter/2, nchan * n_bands, len(rts_baselines)))

#for i in range(n_iter/2):
#    for ch in range(nchan * n_bands):
#        for n in range(len(rts_baselines)):
#            cross_diffs[i,ch,n] = pfb_corrected_vis[2*i,ch,n] - pfb_corrected_vis[(2*i)+1,ch,n]

#cross_rms = np.zeros((nchan * n_bands, len(rts_baselines)))

#for ch in range(nchan * n_bands):
#    for n in range(len(rts_baselines)):
#        cross_rms[ch,n] = np.std(cross_diffs[:,ch,n])

        
