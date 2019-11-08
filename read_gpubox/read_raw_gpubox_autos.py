from matplotlib import pyplot as plt
from load_gpubox_vis import *
import numpy as np
from glob import glob
from spoof_mwac import *
from astropy.io import fits
from fit_cable_reflections import *

n_iter=32
nchan = 32
n_ants = 128
n_bands = 24

# This creates a mapping between the rts baselines (which are indexed according to the metafits 'INPUT' value)

corr_mapping = fill_mapping_matrix()
the_matrix = extract_full_matrix()
rts2index = load_rts_baselines(corr_mapping,the_matrix)
rts_baselines = range(n_ants * (n_ants - 1) / 2)
raw_baselines = [rts2index[i] for i in rts_baselines]
auto2index = load_auto_baselines(corr_mapping,the_matrix)

auto_bl = []


for i in xrange(128):
    if(i==0):
        auto_bl.append(i)
    else:
        auto_bl.append(auto_bl[-1] + i+1)

#raw_autos = zeros((n_iter,nchan,len(auto_bl)),dtype=complex)

mwa_dir = os.getenv('MWA_DIR','/astro/mwaeor/MWA/')

obsid = sys.argv[1]
#meta_file = fits.open(args[1])
#meta_file = fits.open(sys.argv[2])

gpu_files = glob(mwa_dir+'data/'+obsid+'/*gpubox*00.fits')

metafiles = glob(mwa_dir+'data/'+obsid+'/*meta*ppds*fits')
meta_file = fits.open(metafiles[0])
cent_chan = meta_file[0].header['CENTCHAN']
rts_inputs = meta_file[1].data['Input']
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

# need to consider startime in file name(?)    

gpu_files.sort()

band_numbers = [int((g.split('gpubox')[1])[:2])-1 for g in gpu_files]

#gpu_files.reverse()

gpu_files = gpu_files[:n_bands]
                                                                       
raw_autos = zeros((n_iter,nchan*len(gpu_files),len(auto_bl)),dtype=float)


band_0 = cent_chan - 12

#for band,gpubox_file in enumerate(gpu_files):
for band,gpubox_file in zip(band_numbers,gpu_files):

    all_vis = get_raw_vis(gpubox_file,n_iter=n_iter)

    # Correct band ordering
    
    band_number = band + cent_chan - 12
    if(band_number > 128):
        if(band_0 > 128):
            band_index = n_bands -(band+1)
        else:
            band_index = n_bands-(band_number-128)
    else:
        band_index = band

    print band,band_index,gpubox_file
        
    for t in range(n_iter):
        for ch in range(nchan):
            for b in range(len(auto_bl)):
              raw_autos[t,ch+band_index*nchan,b] = all_vis[t][ch][auto2index[b]]  
#              raw_autos[t,ch+band_index*nchan,b] /= pow((float(cable2gains[rts2cables[b]][band_index]) / 64.0),2.0)
                
mean_autos = np.mean(raw_autos,axis=0)
std_autos = np.std(raw_autos,axis=0)

for c in cable_types:
    plt.clf()
    for i in rts_ants:
        if(rts2cables[i] == c):
            if(flags[i]==0):
                plt.plot(mean_autos[:,i])
    plt.savefig('%s_%s.png' % (obsid,c))

#In [197]: band_n = 14

#In [198]: band_medians_14 = np.array([np.median(mean_autos[band_n*nchan:(band_n+1)*nchan,i]) for i in rts_ants if flags[i] == 0])

#In [199]: band_n = 15

#In [200]: band_medians_15 = np.array([np.median(mean_autos[band_n*nchan:(band_n+1)*nchan,i]) for i in rts_ants if flags[i] == 0])

#In [201]: band_n = 16

#In [202]: band_medians_16 = np.array([np.median(mean_autos[band_n*nchan:(band_n+1)*nchan,i]) for i in rts_ants if flags[i] == 0])

#In [203]: plt.plot((band_medians_14 + band_medians_15)/2.0, band_medians_15 / band_medians_14,'ro')
#Out[203]: [<matplotlib.lines.Line2D at 0x7f7f3d538650>]

#In [204]: plt.plot((band_medians_14 + band_medians_15)/2.0, 4.0 * band_medians_15 / band_medians_14,'go')
#Out[204]: [<matplotlib.lines.Line2D at 0x7f7f35703410>]

#In [205]: plt.plot((band_medians_16 + band_medians_15)/2.0, band_medians_16 / band_medians_15,'bo')
#Out[205]: [<matplotlib.lines.Line2D at 0x7f7e029f6cd0>]

#In [206]: plt.clf()

#In [207]: plt.plot((band_medians_14 + band_medians_15)/2.0, band_medians_15 / band_medians_14,'ro',label='Before Gain Jump')
#Out[207]: [<matplotlib.lines.Line2D at 0x7f7e029f63d0>]

#In [208]: plt.plot((band_medians_14 + band_medians_15)/2.0, 4.0 * band_medians_15 / band_medians_14,'go',label='4 x Before Gain Jump')
#Out[208]: [<matplotlib.lines.Line2D at 0x7f7f42d6b050>]

#In [209]: plt.plot((band_medians_16 + band_medians_15)/2.0, band_medians_16 / band_medians_15,'bo',label='Across Gain Jump')

#In [214]: plt.plot((band_medians_16 + band_medians_15)/2.0, band_medians_16 / band_medians_15,'bo',label='Across Gain Jump')

#In [215]: plt.plot((band_medians_14 + 4.0* band_medians_15)/2.0, 4.0 * band_medians_15 / band_medians_14,'go',label='4 x Before Gain Jump')


#In [216]: plt.plot((band_medians_14 + band_medians_15)/2.0, band_medians_15 / band_medians_14,'ro',label='Before Gain Jump')


#In [217]: plt.legend(loc=0)

#In [218]: plt.xlabel('Average of Band Medians')

#In [219]: plt.ylabel('Ratio of Band Medians')


#In [220]: plt.title('XX Autos Before and After Digital Gain Jump for 1061316296')

#In [221]: plt.savefig('XX_Auto_Ratios_1061316296.png')


                                                                                   
        
        
        
