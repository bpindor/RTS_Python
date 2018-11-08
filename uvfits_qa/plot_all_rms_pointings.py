import sys,os
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import max, array

obslist = open(sys.argv[1])

colors = ['r','b','g','c','m']

labelled_pointings = []

for line in obslist:
    obsid = line.split()[0]
    pointing = int(line.split()[1])
    rms_file = '%s_rms.txt' % obsid
    if(os.access(rms_file,os.R_OK)):
        rms = ((open(rms_file)).readline()).split()
        print obsid, (array(rms,dtype=float)[0:len(rms)/2]).max()
        if(pointing not in labelled_pointings):
            plt.plot(rms,colors[pointing+2],label=pointing)
            labelled_pointings.append(pointing)
        else:
            plt.plot(rms,colors[pointing+2])
    else:
        print 'Cannot open file %s' % rms_file

#for e in edges:
#    plt.axvline(e)
#3plt.ylim(0,60)
plt.legend()
plt.xlabel('Fine Channel Number')
plt.ylabel('Visibility Difference RMS')
plt.ylim(0,60)
plt.savefig('all.png')
