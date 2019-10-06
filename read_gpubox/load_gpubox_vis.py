from astropy.io import fits
import os,sys
from numpy import zeros, array, sum


#gpu_file = fits.open(sys.argv[1])

#flagged_channels = [0,1,16,30,31]
#good_channels = [x for x in xrange(32) if x not in flagged_channels]

#auto_bl = []

#for i in xrange(128):
#    if(i==0):
#        auto_bl.append(i)
#    else:
#        auto_bl.append(auto_bl[-1] + i+1)

#na_bl = [i for i in xrange(8256) if i not in auto_bl]

#sum_xx = zeros(27)
#all_xx = zeros((len(gpu_file),27))
#all_bl = zeros((27,8256))
#bl_by_time = []
def get_raw_vis(gpufile,n_iter=32):
    
    gpu_file = fits.open(gpufile)

    if(n_iter > len(gpu_file)):
        n_iter = len(gpu_file)

    all_vis = []

    for i in xrange(1,n_iter+1):
            vis = gpu_file[i].data
            all_vis.append(vis[:,::2] + vis[:,1::2]*1.0j)

    return all_vis
