import RTS_cals
from astropy.io import fits
#from fit_cable_reflections import fit_bandpass_and_reflection, bp_xvals
import os,sys
from optparse import OptionParser,OptionGroup

usage = 'Usage: single_rewrite_cals [directory containing RTS cals] [metafits]'

parser = OptionParser(usage=usage)

parser.add_option('--n_bands',dest='n_bands',type='int',default=24,help='Number of Coarse Bands Present')

(options, args) = parser.parse_args()

cal_dir = args[0]
metafile = args[1]

n_bands = options.n_bands

raw_cal = RTS_cals.rts_cal(n_bands=n_bands)
fit_cal = RTS_cals.rts_cal(n_bands=n_bands)
copy_cal = RTS_cals.rts_cal(n_bands=n_bands)
copy_fit = RTS_cals.rts_cal(n_bands=n_bands)
raw_cal.load_all_BP_jones(path=cal_dir, raw=True,n_bands=n_bands)
fit_cal.load_all_BP_jones(path=cal_dir, raw=False,n_bands=n_bands)
raw_cal.load_all_DI_jones(path=cal_dir,n_bands=n_bands)
fit_cal.load_all_DI_jones(path=cal_dir,n_bands=n_bands)
raw_cal.form_single_jones()
fit_cal.form_single_jones()

meta_file = fits.open(metafile)
cent_chan = meta_file[0].header['CENTCHAN']

for i,a in enumerate(raw_cal.antennas):
    if(raw_cal.antennas[i].BP_jones[0] is not None):
        single_jones = raw_cal.all_single_jones(antenna=i,cent_chan=cent_chan)
        copy_cal.load_cals_from_single_jones(single_jones,antenna=i,cent_chan=cent_chan)
        single_jones_fit = fit_cal.all_single_jones(antenna=i,cent_chan=cent_chan)
        copy_fit.load_cals_from_single_jones(single_jones_fit,antenna=i,cent_chan=cent_chan)
        
#        single_jones_fit = fit_bandpass_and_reflection(single_jones)
#        fit_cal.load_cals_from_single_jones(single_jones_fit,antenna=i,cent_chan=cent_chan)

RTS_cals.write_BP_files(raw_cal,fit_cal,filename='fit_test')
RTS_cals.write_BP_files(copy_cal,copy_fit,filename='copy_test')
RTS_cals.write_DI_files(copy_cal,filename='copy_test')

for i in range(n_bands):
    cmd = 'cp Bandpass_copy_test%03d.dat BandpassCalibration_node%03d.dat' % (i+1, i+1)
    print cmd
    os.system(cmd)
    cmd = 'cp DI_Jones_copy_test%03d.dat DI_JonesMatrices_node%03d.dat' % (i+1, i+1)
    print cmd
    os.system(cmd)

#copy2_cal = RTS_cals.rts_cal()
#copy2_cal.load_all_BP_jones(path='./', raw=False)
#copy2_cal.load_all_DI_jones(path='./')
#copy2_cal.form_single_jones()

#single_jones_copy = copy2_cal.all_single_jones(antenna=0,cent_chan=cent_chan)
    
    
