#! /usr/bin/env python

#####

def generate_galaxy(infile):

    n_obs = sum(1 for line in open(infile))
 
    out_file = open('sReflag_mwaf_files.sh','w+')

    out_file.write('#!/bin/bash -l\n')
    out_file.write('#SBATCH --nodes=1\n')
    out_file.write('#SBATCH --job-name=\"reflag_mwaf\"\n')
    out_file.write('#SBATCH -o reflag_mwaf-%A-%a.out\n')
    out_file.write('#SBATCH --array=1-%s\n' % n_obs)
    out_file.write('#SBATCH --ntasks-per-node=1\n')
    out_file.write('#SBATCH --time=00:90:00\n')
    out_file.write('#SBATCH --partition=gpuq\n')
    out_file.write('#SBATCH --account=mwaeor\n')
    out_file.write('#SBATCH --export=NONE\n')

    # out_file.write('source /scratch2/mwaeor/bpindor/MWA_Python/bin/activate\n')
    out_file.write('cd $SLURM_SUBMIT_DIR\n')
    out_file.write('aprun reflag_mwaf_files.py %s \n' % infile)

    out_file.close()

def generate_ozstar(infile):

    n_obs = sum(1 for line in open(infile))
 
    out_file = open('sReflag_mwaf_files.sh','w+')

    out_file.write('#!/bin/bash -l\n')
    out_file.write('#SBATCH --nodes=1\n')
    out_file.write('#SBATCH --job-name=\"reflag_mwaf\"\n')
    out_file.write('#SBATCH -o reflag_mwaf-%A-%a.out\n')
    out_file.write('#SBATCH --array=1-%s\n' % n_obs)
    out_file.write('#SBATCH --ntasks-per-node=1\n')
    out_file.write('#SBATCH --time=00:30:00\n')
    out_file.write('#SBATCH --partition=skylake\n')
    out_file.write('#SBATCH --account=oz048\n')
    out_file.write('#SBATCH --mem=5000\n')
    out_file.write('#SBATCH --cpus-per-task=2\n')
    out_file.write('#SBATCH --export=NONE\n')

    
    out_file.write('cd $SLURM_SUBMIT_DIR\n')
    out_file.write('srun --export=ALL --mem=5000  reflag_mwaf_files.py %s \n' % infile)

    out_file.close()
    

#####

import sys,os, glob
from optparse import OptionParser,OptionGroup
from astropy.io import fits

usage = 'generate_sReflag_mwaf.py [text file of obsIDs]'

parser = OptionParser(usage=usage)

(options, args) = parser.parse_args()

infile = args[0]

mwa_dir = os.getenv('MWA_DIR','/scratch/partner1019/MWA/')

if(mwa_dir == '/astro/mwaeor/MWA/'):
    generate_galaxy(infile)
if(mwa_dir == '/fred/oz048/MWA/'):
    generate_ozstar(infile)
