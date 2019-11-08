#! /usr/bin/env python
"""
generate_obsIDs2mwacRTS_multi.py
Takes a list of obsIDs and generates a qsub script for downloading the data from NGAS (if requested), running ngas2mwaRTS.py and another qsub script for running the RTS  

09/08/13 - Added different functions for generating scripts on gstar and fornax given their quite different I/O structure
"""

##########

def generate_gstar(infile,basename,n_subbands,rts_templates,download,array,options):

    mwa_dir = os.getenv('MWA_DIR','/scratch/astronomy556/MWA/')
    cwd = os.getcwd()

    n_obs = sum(1 for line in open(infile))

    out_file = open('NGAS2mwacRTS_multi.sh','w+')
    
    out_file.write('#!/bin/bash\n')
    for line in in_file:
        obs_id = (line.split())[0]
        if(download):
            out_file.write('cd ' + mwa_dir + 'data\n')
            out_file.write('change_db.py curtin\n')
            out_file.write('obsresolve.py -r ngas01.ivec.org -s ngas01.ivec.org -o %s\n' % obs_id)
        data_dir = mwa_dir + 'data/%s' % obs_id    
        out_file.write('cd '+ data_dir +'\n')
        out_file.write('change_db.py curtin\n') 

        if(array == '32T'):
            out_file.write('ngas2uvfitssubbandsRTS.py %s %s %s %s %s\n' % (basename, data_dir, obs_id, n_subbands, rts_templates))
        else:
            out_file.write('ngas2mwacRTS.py %s %s %s %s %s\n' % (basename, data_dir, obs_id, n_subbands, rts_templates))

    # Command to write RTS qsub script

    out_file.write('cd %s\n' %cwd)
    out_file.write('generate_mwac_qRTS.py %s %s %s %s \n' % (infile, basename, n_subbands, rts_templates))

    out_file.close()
    


def generate_fornax(infile,basename,n_subbands,rts_templates,download,array,options):

    mwa_dir = os.getenv('MWA_DIR','/scratch/astronomy556/MWA/')
    cwd = os.getcwd()

    n_obs = sum(1 for line in open(infile))

    if options.do_erase is False:
        do_erase = ''
    else:
        do_erase = '--erase_GPUfiles'

    if(options.channel_flags != mwa_dir + '/RTS/utils/flags/flagged_channels_default.txt'):
        flags_string = '--channel_flags=' + options.channel_flags
    else:
        flags_string = ''

    if(options.dynamic is True):
        dynamic_string = '--dynamic_sourcelist=%d' % options.dynamic 
    else:
        dynamic_string = ''
    if(options.tile_flags is not None):
        tile_flags_string = '--tile_flags=' + options.tile_flags
    else:
        tile_flags_string = ''

    if(options.process_time):
        process_string = '--process_time=%d' % options.process_time
    else:
        process_string = ''

    out_file = open('qNGAS2mwacRTS_multi.sh','w+')

    out_file.write('#!/bin/bash\n')
    out_file.write('#PBS -l select=1:ncpus=1:mem=8gb\n')
    out_file.write('#PBS -l walltime=00:%d:00\n' % (n_obs * 30))
    out_file.write('#PBS -m e\n')
    out_file.write('#PBS -q copyq\n')
    out_file.write('#PBS -W group_list=partner1019\n')
    out_file.write('source ' + mwa_dir + 'bin/activate \n')


    for line in in_file:
        obs_id = (line.split())[0]
        if(download):
            out_file.write('cd ' + mwa_dir + 'data\n')
            #out_file.write('change_db.py curtin\n')
            out_file.write('obsresolve.py -r ngas01.ivec.org -s ngas01.ivec.org -o %s\n' % obs_id)
            data_dir = mwa_dir + 'data/%s' % obs_id    
            out_file.write('cd '+ data_dir +'\n')

        
            if(array == '32T'):
                out_file.write('ngas2uvfitssubbandsRTS.py %s %s %s %s %s\n' % (basename, data_dir, obs_id, n_subbands, rts_templates))
            else:
                out_file.write('ngas2mwacRTS.py %s %s %s %s %s\n' % (basename, data_dir, obs_id, n_subbands, rts_templates))

    # Command to write RTS qsub script

    out_file.write('cd $PBS_O_WORKDIR\n')
    if(options.auto_gen):
        out_file.write('generate_mwac_qRTS_auto.py %s %s %s %s --auto --chunk_number=%d %s %s %s %s %s\n' % (infile, basename, n_subbands, rts_templates, options.chunk_number,do_erase,flag_string,dynamic_string,tile_flags_string,process_string))
    else:
        if options.use_metafits:
            out_file.write('generate_mwac_qRTS.py %s %s %s %s \n' % (infile, basename, n_subbands, rts_templates))
        else:
            out_file.write('generate_mwac_qRTS.py %s %s %s %s \n' % (infile, basename, n_subbands, rts_templates))

    out_file.write('deactivate\n')

    out_file.close()

def generate_galaxy(infile,basename,rts_templates,options):

    mwa_dir = os.getenv('MWA_DIR','/astro/mwaeor/MWA/')
    cwd = os.getcwd()

    n_obs = sum(1 for line in open(infile))

    if options.do_erase is False:
        do_erase = ''
    else:
        do_erase = '--erase_GPUfiles'

    if(options.channel_flags != mwa_dir + '/RTS/utils/flags/flagged_channels_default.txt'):
        flags_string = '--channel_flags=' + options.channel_flags
    else:
        flags_string = ''

    if(options.dynamic):
        dynamic_string = '--dynamic_sourcelist=%d' % options.dynamic
    else:
        dynamic_string = ''

    if(options.sourcelist):
        sourcelist_string =  '--sourcelist=%s' % options.sourcelist
    else:
        sourcelist_string =  '--sourcelist=/astro/mwaeor/bpindor/PUMA/srclists/srclist_puma-v2_complete.txt'

    if(options.tile_flags is not None):
        tile_flags_string = '--tile_flags=' + options.tile_flags
    else:
        tile_flags_string = ''

    if options.dev_rts:
        dev_rts_string = '--dev_rts'
    else:
        dev_rts_string= ''

    if options.tag_uvfits:
        tag_uvfits_string = '--tag_uvfits'
    else:
        tag_uvfits_string= ''

    if options.tag_logs:
        tag_logs_string = '--tag_logs'
    else:
        tag_logs_string= ''

    if(options.process_time):
        process_string = '--process_time=%d' % options.process_time
    else:
        process_string = ''        


    out_file = open('sGetNGASData_pipeline_%s.sh' % options.chunk_number,'w+')
    
    out_file.write('#!/bin/bash -l\n')
    out_file.write('#SBATCH --job-name=\"getNGASdata\"\n')
    out_file.write('#SBATCH -o getNGASdata-%A.out \n')
    out_file.write('#SBATCH --ntasks=1\n')
    out_file.write('#SBATCH --ntasks-per-node=1\n')
    out_file.write('#SBATCH --time=00:%d:00\n' % (n_obs * 30))
    out_file.write('#SBATCH --account=mwaeor\n')
    out_file.write('#SBATCH --export=NONE\n')
    if(options.rts_only):
        out_file.write('#SBATCH --partition=gpuq\n') 
    else:
        out_file.write('#SBATCH --clusters=zeus\n')
        out_file.write('#SBATCH --partition=copyq\n')

        out_file.write('module load pyephem\n')
        out_file.write('module load setuptools\n')
        out_file.write('cd ' + mwa_dir + 'data\n')
        for line in in_file:
            obs_id = (line.split())[0]
            # out_file.write('cd ' + mwa_dir + 'data\n')
            # out_file.write('obsdownload.py -o %s \n' % obs_id)
            # out_file.write('cd ' + mwa_dir + 'data/%s \n' % obs_id)
            # out_file.write('make_metafits.py -g %s \n' % obs_id)

            if (options.no_gpubox == False):
                out_file.write('obsdownload.py -o %s --chstart=1 --chcount=24\n' % obs_id)
            out_file.write('obsdownload.py -o %s -f\n' % obs_id)

            out_file.write('obsdownload.py -m -o %s \n' % obs_id)

            
    out_file.write('cd $SLURM_SUBMIT_DIR\n')
    if (options.no_gpubox == False):
        out_file.write('list_gpubox_files.py %s\n' % infile)
        out_file.write('./autoProcess_pipeline_%d.sh\n' % options.chunk_number)

    out_file.close()

def generate_ozstar(infile,basename,rts_templates,options):

    n_obs = sum(1 for line in open(infile))

    out_file = open('sGetNGASData_pipeline_%s.sh' % options.chunk_number,'w+')

    out_file.write('#!/bin/bash -l\n')
    out_file.write('#SBATCH --job-name=\"getNGASdata\"\n')
    out_file.write('#SBATCH -o getNGASdata-%A.out \n')
    out_file.write('#SBATCH --ntasks=1\n')
    out_file.write('#SBATCH --ntasks-per-node=1\n')
    out_file.write('#SBATCH --time=00:%d:00\n' % (n_obs * 30))
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

usage = 'Usage: generate_getNGASData.py [text file of obsIDs] [basename] [# of subbands] [rts_templates]'

parser = OptionParser(usage=usage)

parser.add_option('--auto',action="store_true",dest="auto_gen", help="Generate script to be run automatically from within qRTS_multi_wrapper")
parser.add_option('--erase_GPUfiles',action='store_true',dest='do_erase',default=False,help='Erase raw GPU files after processing [default=%default]')
parser.add_option('--use_metafits',action='store_true',dest='use_metafits',default=True,help='Use metafits file to gather metadata [default=%default]')
parser.add_option('--chunk_number',dest='chunk_number',type='int',default=0,
                  help='Chunk number for automatic processing')
parser.add_option('--channel_flags',dest="channel_flags",default=mwa_dir + '/RTS/utils/flags/flagged_channels_default.txt',
                      help="File containing flagged channels",metavar="CHANNELFLAGS") 
parser.add_option('--dynamic_sourcelist',dest='dynamic',type='int',default=0,help='Use Dynamically Generated Sourcelist for Calibration [number of sources] (only applies to first entry in template file)')    
parser.add_option('--sourcelist',dest='sourcelist',type='string',default='',help='Specify base catalog for dynamic sourcelist (only applies to first entry in template file)')
parser.add_option('--tile_flags',dest="tile_flags",default=None,
                      help="Comma separated list of additional tiles to flag. Use tile names e.g. 8th tile on Rx10 is Tile108, so --tile_flags=108. Default is none.")
parser.add_option('--dev_rts',dest="dev_rts",default=False,action='store_true',
                  help="Use development version of the RTS")
parser.add_option('--tag_uvfits',dest="tag_uvfits",default=False,action='store_true',
                  help="Tag and copy uvfits to data subdirectory")
parser.add_option('--tag_logs',dest="tag_logs",default=False,action='store_true',
                  help="Tag and copy uvfits to data subdirectory")
parser.add_option('--process_time',dest="process_time",type='int',default=0,help='Processing time per RTS run per obsid. Default is 20 minutes.')
parser.add_option('--rts_only',dest="rts_only",default=False,action='store_true',help="Skip downloading and RTS setup steps")
parser.add_option('--no_gpubox',dest="no_gpubox",default=False,action='store_true',help="Only download flag and metafits files")



(options, args) = parser.parse_args()

infile = args[0]
basename = args[1]
rts_templates = args[2]



try:
    in_file = open(infile)

except IOError, err:
    'Cannot open input file %s\n',str(infile)

mwa_dir = os.getenv('MWA_DIR','/scratch/partner1019/MWA/')

if(mwa_dir == '/lustre/projects/p048_astro/MWA/'):
    generate_gstar(infile,basename,n_subbands,rts_templates,download,array,options)
if(mwa_dir == '/scratch/partner1019/MWA/'):
    generate_fornax(infile,basename,n_subbands,rts_templates,download,array,options)
if(mwa_dir == '/astro/mwaeor/MWA/'):
    generate_galaxy(infile,basename,rts_templates,options)

if(mwa_dir == '/fred/oz048/MWA/'):
    generate_ozstar(infile,basename,rts_templates,options)
        


    
    
    



    
        
        
