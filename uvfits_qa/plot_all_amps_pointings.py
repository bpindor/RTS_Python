import sys,os
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import max, array

obslist = open(sys.argv[1])

colors = ['r','b','g','c','m']

labelled_pointings = []

pol = sys.argv[2]

for line in obslist:
    obsid = line.split()[0]
    pointing = int(line.split()[1])
    amps_file = '%s_%s_amps.txt' % (obsid,pol)
    if(os.access(amps_file,os.R_OK)):
        amps = ((open(amps_file)).readline()).split()
        print obsid, (array(amps,dtype=float)[0:len(amps)/2]).max()
        if(pointing not in labelled_pointings):
            plt.plot(amps,colors[pointing+2],label=pointing)
            labelled_pointings.append(pointing)
        else:
            plt.plot(amps,colors[pointing+2])
    else:
        print 'Cannot open file %s' % amps_file

#for e in edges:
#    plt.axvline(e)
#3plt.ylim(0,60)
#plt.legend()
plt.xlabel('Fine Channel Number')
plt.ylabel('Amps (XX)')
plt.ylim((0,60.0))
#plt.legend()
plt.savefig('all_%s_pointings.png' % pol)
