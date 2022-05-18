import RTS_cals
from matplotlib import pyplot as plt

obs_list = ['1125953008']
mwa_dir = '/home/bpindor/code/MWA/'
subdir = 'EOR1_2015_baseline'
n_bands = 24
metafile = '/home/bpindor/code/MWA//data/1125953008/1125953008_metafits_ppds.fits'


for obs in obs_list:

    data_path = f'{mwa_dir}data/{obs}/{subdir}/'

    raw_cal = RTS_cals.rts_cal(n_bands=n_bands)

    raw_cal.load_all_cals(path=data_path,metafile=metafile)

# Plot XX Amps of a given antenna
    
antenna_n = 0

plt.plot(abs(raw_cal.antennas[antenna_n].full_band_jones[:,0,0]))
plt.savefig('XX_Amps_Antenna_%d.png' % antenna_n)
    

    
