#! /usr/bin/env python

""" 
generate_getGPUBOX_Data.py
Takes a list of obsids and generates a script which launches the RTS processing pipeline. Eventually should wrap downloading option via mantra ray client
"""

def generate_ozstar(infile,options):

    n_obs = sum(1 for line in open(infile))

    out_file = open('sgetGPUBOXData_pipeline_%s.sh' % options.chunk_number,'w+')

    out_file.write('#!/bin/bash -l\n')
    out_file.write('#SBATCH --job-name=\"getGPUBOXdata\"\n')
    out_file.write('#SBATCH -o getGPUBOXdata-%A.out \n')
    out_file.write('#SBATCH --ntasks=1\n')
    out_file.write('#SBATCH --ntasks-per-node=1\n')
    out_file.write('#SBATCH --time=00:%d:00\n' % (n_obs * 5))
    out_file.write('#SBATCH --account=oz048\n')
    out_file.write('#SBATCH --partition=skylake\n')
    out_file.write('#SBATCH --export=NONE\n')

    
    if (options.no_gpubox == False):
        out_file.write('srun --export=ALL list_gpubox_files.py %s\n' % infile)
        out_file.write('./autoProcess_pipeline_%d.sh\n' % options.chunk_number)

    out_file.close()

    
    

##########

import sys,os, glob
from optparse import OptionParser,OptionGroup

mwa_dir = os.getenv('MWA_DIR','/astro/mwaeor/MWA/')

usage = 'Usage: generate_getGPUBOX_Data.py [text file of obsIDs]'

parser = OptionParser(usage=usage)

parser.add_option('--chunk_number',dest='chunk_number',type='int',default=0,
                  help='Chunk number for automatic processing')
parser.add_option('--no_gpubox',dest="no_gpubox",default=False,action='store_true',help="Only download flag and metafits files")

(options, args) = parser.parse_args()

infile = args[0]

try:
    in_file = open(infile)

except IOError, err:
    'Cannot open input file %s\n',str(infile)

mwa_dir = os.getenv('MWA_DIR','/scratch/partner1019/MWA/')

if(mwa_dir == '/astro/mwaeor/MWA/'):
    generate_galaxy(infile,options)

if(mwa_dir == '/fred/oz048/MWA/'):
    generate_ozstar(infile,options)
