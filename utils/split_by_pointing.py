import sys,os
import glob
import astropy.io.fits as fits
from optparse import OptionParser,OptionGroup

usage = 'Usage: split_by_pointing.py [text file of obsIDs]'

parser = OptionParser(usage=usage)

(options, args) = parser.parse_args()

infile = args[0]

try:
    in_file = open(infile)

except IOError, err:
    'Cannot open input file %s\n',str(infile)

mwa_dir = os.getenv('MWA_DIR','/scratch/astronomy556/MWA/')

grid_files = {}
grid_list = []
#file_list = []
out_file = open('all_pointings.txt','w+')

for line in in_file:
    obsid = (line.split())[0]
    data_dir = mwa_dir + 'data/' + obsid
    metafile = (glob.glob(data_dir + '/*meta*fits'))[0]
    meta_file = fits.open(metafile)
    gridnum = int((meta_file[0].header)['GRIDNUM'])
    if (gridnum not in grid_list):
        grid_list.append(gridnum)
        grid_files[gridnum] = open('pointing_%d.dat' % gridnum,'w+')
    grid_files[gridnum].write('%s\n' % obsid)
    out_file.write('%s %d\n' % (obsid, gridnum))
    #print (meta_file[0].header)['GRIDNUM']

out_file.close()                  
for num in grid_files.keys():
    grid_files[num].close()
