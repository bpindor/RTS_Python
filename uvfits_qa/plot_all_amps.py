import sys,os
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import max, array

obslist = open(sys.argv[1])

for line in obslist:
    obsid = line.split()[0]
    amps_file = '%s_amps.txt' % obsid
    if(os.access(amps_file,os.R_OK)):
        amps = ((open(amps_file)).readline()).split()
        amps_min = array([r for r in amps if float(r) > 0.0],dtype=float).min()
        print obsid, (array(amps,dtype=float)[0:len(amps)/2]).max(),amps_min
        plt.plot(array([float(r) for r in amps]))
    else:
        print 'Cannot open file %s' % amps_file

#for e in edges:
#    plt.axvline(e)
#plt.ylim(50,250)
#plt.ylim(10,40)
plt.savefig('all_amps.png')
