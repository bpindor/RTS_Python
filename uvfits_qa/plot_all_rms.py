import sys,os
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import max, array

obslist = open(sys.argv[1])

plot_min = 1.0e6
plot_max = 0.0

for line in obslist:
    obsid = line.split()[0]
    rms_file = '%s_rms.txt' % obsid
    if(os.access(rms_file,os.R_OK)):
        rms = ((open(rms_file)).readline()).split()
        rms_min = array([r for r in rms if float(r) > 0.0],dtype=float).min()
        rms_max = array([r for r in rms if float(r) > 0.0],dtype=float).max()
        if(rms_min < plot_min):
            plot_min = rms_min
        if(rms_max > plot_max):
            plot_max = rms_max
        print( obsid, (array(rms,dtype=float)[0:int(len(rms)/2)]).max(),rms_min)
        plt.plot(array([float(r) for r in rms]))
    else:
        print('Cannot open file %s' % rms_file)

#for e in edges:
#    plt.axvline(e)
#plt.ylim(50,250)
plt.ylim(plot_min/1.2,plot_max*1.2)
plt.savefig('all.png')
