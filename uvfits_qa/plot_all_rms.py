import sys,os
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import max, array

obslist = open(sys.argv[1])

for line in obslist:
    obsid = line.split()[0]
    rms_file = '%s_rms.txt' % obsid
    if(os.access(rms_file,os.R_OK)):
        rms = ((open(rms_file)).readline()).split()
        rms_min = array([r for r in rms if float(r) > 0.0],dtype=float).min()
        print obsid, (array(rms,dtype=float)[0:len(rms)/2]).max(),rms_min
        plt.plot(array([float(r) for r in rms]))
    else:
        print 'Cannot open file %s' % rms_file

#for e in edges:
#    plt.axvline(e)
#plt.ylim(50,250)
plt.ylim(10,40)
plt.savefig('all.png')
