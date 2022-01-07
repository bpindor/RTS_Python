#! /usr/bin/env python
"""
generate_sRTS_auto.py
Generates an sbatch script to process a series of list of obsIDs. This script is run automatically at the end of ??.sh
04/09/13: Writes a two part 

"""
import string 

def generate_galaxy(infile,basename,n_subbands,rts_templates,array,options):
    try:
        in_file = open(infile)

    except IOError, err:
        'Cannot open input file %s\n',str(infile)

    n_obs = sum(1 for line in open(infile))
    n_templates = sum(1 for line in open(rts_templates))


    if (options.auto_gen):
        rts_file = open('sRTS_auto_wrapper_%d.sh' % options.chunk_number,'w+')
        rts_file.write('#!/bin/bash -l\n')
        rts_file.write('#SBATCH --job-name=\"RTS\"\n')
        rts_file.write('#SBATCH -o RTS-%A.out\n')
        rts_file.write('#SBATCH --nodes=%d\n' % (int(n_subbands) + 1))
        rts_file.write('#SBATCH --ntasks-per-node=1\n')
        rts_file.write('#SBATCH --time=00:%d:00\n' % (n_obs * n_templates * options.process_time))
        rts_file.write('#SBATCH --partition=gpuq\n')
        rts_file.write('#SBATCH --account=mwaeor\n')
        rts_file.write('#SBATCH --export=NONE\n')
#        rts_file.write('#SBATCH --mem=32000\n')
#        rts_file.write('source /scratch2/mwaeor/bpindor/MWA_Python/bin/activate \n')
        rts_file.write('cd $SLURM_SUBMIT_DIR\n')
        rts_file.write('./sRTS_auto_inner_%d.sh\n' % options.chunk_number)

        rts_file.close() 

        rts_file = open('sRTS_auto_inner_%d.sh' % options.chunk_number,'w+')

        if(options.channel_flags != mwa_dir + '/RTS/utils/flags/flagged_channels_default.txt'):
            flags_string = '--channel_flags=' + options.channel_flags
        else:
            flags_string = ''

        if(options.tile_flags is not None):
            tile_flags_string = '--tile_flags=' + options.tile_flags
        else:
            tile_flags_string = ''

        if options.dev_rts:
            rts_string = '/group/mwaeor/bpindor/RTS/bin/rts_gpu'
        else:
            rts_string= 'rts_gpu'

        if(options.sourcelist):
           sourcelist =  options.sourcelist
           sourcelist_basename=string.split(sourcelist,'.txt')[0]
           sourcelist_basename=string.split(sourcelist_basename,'/')[-1]
        else:
           sourcelist_basename='srclist_puma-v2_complete'

        for line in in_file:
            obs_id = (line.split())[0]
            data_dir = mwa_dir + 'data/%s' % obs_id    
            rts_file.write('cd '+ data_dir +'\n')

            try:
                template_list_file = open(rts_templates)
            except IOError, err:
                'Cannot open list of RTS template files %s\n' % rts_templates
            #check if old-style metafits file exists:
            metafits_filename='%s.metafits' % obs_id
            metafits_fullpath='%sdata/%s/%s.metafits' % (mwa_dir,obs_id,obs_id)
            if os.path.isfile(metafits_fullpath):
               pass
               print 'using old-style metafits file'
            else:
               metafits_filename='%s_metafits_ppds.fits' % obs_id
            rts_file_index = 0
            rts_file.write('generate_RTS_in_mwac.py %s %s %s %s --templates=%s --header=%s %s %s\n' % (data_dir, basename, n_subbands,array,rts_templates,metafits_filename,flags_string,tile_flags_string))
            if(options.dynamic):
                #rts_file.write('sed -i s,/scratch2/mwaeor/bpindor/pp_selfcal/uniq_300conv_eor0.txt,%s/srclist_puma-v2_complete_%s_patch%d.txt,g %s_rts_0.in\n' % (data_dir,obs_id,options.dynamic,basename))
                rts_file.write('sed -i s,/scratch2/mwaeor/bpindor/pp_selfcal/uniq_300conv_eor0.txt,%s/%s_%s_patch%d.txt,g %s_rts_0.in\n' % (data_dir,sourcelist_basename,obs_id,options.dynamic,basename))
            # Account for missing correlator GPU files
            t00_files = glob.glob(data_dir + '/*gpubox*00.fits') 
            n_presentBands = len(t00_files)

            for line in template_list_file:
                rts_file.write('aprun -n %d -N 1 %s %s_rts_%d.in > aprun.${SLURM_JOB_ID}_%d.log\n' % (int(n_presentBands)+1, rts_string, basename, rts_file_index,rts_file_index))
                rts_file_index += 1
            template_list_file.close()

            # Clean up
            if options.do_erase:
                rts_file.write('rm *gpubox*.fits\n')

        # Ends loop over obsids
        
        if(options.tag_uvfits):
            rts_file.write('tag_uvfits.py %s %s\n' % (infile, basename))
        if(options.tag_logs):
            rts_file.write('tag_RTS_logs.py %s ${SLURM_JOB_ID} %s\n' % (infile, basename))
        rts_file.close()

def generate_ozstar(infile,basename,n_subbands,rts_templates,array,options):
    try:
        in_file = open(infile)

    except IOError, err:
        'Cannot open input file %s\n',str(infile)

    n_obs = sum(1 for line in open(infile))
    n_templates = sum(1 for line in open(rts_templates))

    if (options.auto_gen):

        rts_file = open('sRTS_auto_wrapper_%d.sh' % options.chunk_number,'w+')
        rts_file.write('#!/bin/bash -l\n')
        rts_file.write('#SBATCH --job-name=\"RTS\"\n')
        rts_file.write('#SBATCH -o RTS-%A.out\n')
        rts_file.write('#SBATCH --nodes=%d\n' % (int(n_subbands)/options.gpus_per_node + 1))
        rts_file.write('#SBATCH --ntasks-per-node=%d\n' % options.gpus_per_node)
        rts_file.write('#SBATCH --time=00:%d:00\n' % (n_obs * n_templates * options.process_time))
        rts_file.write('#SBATCH --partition=skylake-gpu\n')
        rts_file.write('#SBATCH --account=oz048\n')
        rts_file.write('#SBATCH --export=NONE\n')
        rts_file.write('#SBATCH --mem=%d\n' % options.mem)
        rts_file.write('#SBATCH --gres=gpu:%d\n' % options.gpus_per_node)
#        rts_file.write('#SBATCH --cpus-per-task=2\n')
        
        rts_file.write('cd $SLURM_SUBMIT_DIR\n')
        rts_file.write('./sRTS_auto_inner_%d.sh\n' % options.chunk_number)

        rts_file.close() 
        

        rts_file = open('sRTS_auto_inner_%d.sh' % options.chunk_number,'w+')

        if(options.channel_flags != mwa_dir + '/RTS/utils/flags/flagged_channels_default.txt'):
            flags_string = '--channel_flags=' + options.channel_flags
        else:
            flags_string = ''

        if(options.tile_flags is not None):
            tile_flags_string = '--tile_flags=' + options.tile_flags
        else:
            tile_flags_string = ''

        if(options.sourcelist):
           sourcelist =  options.sourcelist
           sourcelist_basename=string.split(sourcelist,'.txt')[0]
           sourcelist_basename=string.split(sourcelist_basename,'/')[-1]
        else:
           sourcelist_basename='srclist_puma-v2_complete'
            
        if options.dev_rts:
            rts_string = '/fred/oz048/bpindor/mwa-RTS/bin/rts_gpu'
            #rts_string = '/fred/oz048/bpindor/RTS_alt/bin/rts_gpu'
        else:
            rts_string= '/fred/oz048/MWA/CODE/mwa-RTS/bin/rts_gpu'

           
        for line in in_file:
            obs_id = (line.split())[0]
            data_dir = mwa_dir + 'data/%s' % obs_id

            try:
                template_list_file = open(rts_templates)
            except IOError, err:
                'Cannot open list of RTS template files %s\n' % rts_templates

            #check if old-style metafits file exists:
            metafits_filename='%s.metafits' % obs_id
            metafits_fullpath='%sdata/%s/%s.metafits' % (mwa_dir,obs_id,obs_id)
            if os.path.isfile(metafits_fullpath):
               pass
               print 'using old-style metafits file'
            else:
               metafits_filename='%s_metafits_ppds.fits' % obs_id
               metafits_fullpath='%sdata/%s/%s' %(mwa_dir,obs_id,metafits_filename)
    
            rts_file_index = 0            
            rts_file.write('cd '+ data_dir + '\n')
            rts_file.write('mkdir -p %s\n' % basename)
            rts_file.write('cp gpufiles_list.dat %s\n' % basename)
            rts_file.write('cd %s\n' % basename)
            rts_file.write('generate_RTS_in_auto.py %s %s %s %s --templates=%s --header=%s %s %s\n' % (data_dir, basename, n_subbands,array,rts_templates,metafits_fullpath,flags_string,tile_flags_string))
            if(options.dynamic):
                rts_file.write('sed -i s,/group/mwaeor/bpindor/pp_selfcal/uniq_300conv_eor0.txt,%s/%s_%s_patch%d.txt,g %s_rts_0.in\n' % (data_dir,sourcelist_basename,obs_id,options.dynamic,basename))
                #rts_file.write('sed -i s,Replace_name_of_sourcecatalogFile,%s/%s_%s_patch%d.txt,g %s_rts_0.in\n' % (data_dir,sourcelist_basename,obs_id,options.dynamic,basename))
            for line in template_list_file:
                rts_file.write('srun --export=ALL --mem=%d --ntasks=%d  --nodes=%d --gres=gpu:%d  --ntasks-per-node=%d %s %s_rts_%d.in > srun.${SLURM_JOB_ID}_%d.log\n' % (options.mem, int(n_subbands) + 1,int(n_subbands)/options.gpus_per_node + 1, options.gpus_per_node,options.gpus_per_node, rts_string,basename,rts_file_index,rts_file_index))
                rts_file_index += 1
            template_list_file.close()

        rts_file.close()
###

import sys,os, glob
from optparse import OptionParser,OptionGroup

usage = 'Usage: generate_sRTS_auto.py [text file of obsIDs] [basename] [# of subbands] [rts_templates]'

mwa_dir = os.getenv('MWA_DIR','/scratch/partner1019/MWA/')

parser = OptionParser(usage=usage)

parser.add_option('--auto',action="store_true",dest="auto_gen", help="Generate script to be run automatically from within qRTS_multi_wrapper")
parser.add_option('--no_upload',action='store_true',dest='no_upload',default=False,help='Do not upload cal solution to eorlive webpage [default=%default]')
parser.add_option('--erase_GPUfiles',action='store_true',dest='do_erase',default=False,help='Erase raw GPU files after processing [default=%default]')
parser.add_option('--chunk_number',dest='chunk_number',type='int',default=0,
                  help='Chunk number for automatic processing')
parser.add_option('--channel_flags',dest="channel_flags",default=mwa_dir + '/RTS/utils/flags/flagged_channels_default.txt',
                      help="File containing flagged channels",metavar="CHANNELFLAGS")
parser.add_option('--dynamic_sourcelist',dest='dynamic',type='int',default=0,help='Use Dynamically Generated Sourcelist for Calibration [number of sources]')
parser.add_option('--sourcelist',dest='sourcelist',type='string',default='',help='Specify the base catalog to be used to generate the dynamic sourcelist (only applies to first entry in template file) default is /scratch2/mwaeor/bpindor/PUMA/srclists/srclist_puma-v2_complete.txt')
parser.add_option('--tile_flags',dest="tile_flags",default=None,
                      help="Comma separated list of additional tiles to flag. Use tile names e.g. 8th tile on Rx10 is Tile108, so --tile_flags=108. Default is none.")
parser.add_option('--dev_rts',dest="dev_rts",default=False,action='store_true',
                  help="Use development version of the RTS")
parser.add_option('--tag_uvfits',dest="tag_uvfits",default=False,action='store_true',
                  help="Tag and copy uvfits to data subdirectory")
parser.add_option('--tag_logs',dest="tag_logs",default=False,action='store_true',
                  help="Tag and copy uvfits to data subdirectory")
parser.add_option('--process_time',dest="process_time",type='int',default=20,help='Processing time per RTS run per obsid. Default is 20 minutes.')
parser.add_option('--gpus_per_node',dest="gpus_per_node",type='int',default=1,help='Number of GPUS to use per node')
parser.add_option('--mem',dest="mem",type='int',default=5000,help='RAM requested for RTS per node. In MB. Default 5000.')

(options, args) = parser.parse_args()

infile = args[0]
basename = args[1]
n_subbands = int(args[2])
rts_templates = args[3]
array = '128T' # No 32T option for this script

if(mwa_dir == '/astro/mwaeor/MWA/'):
    generate_galaxy(infile,basename,n_subbands,rts_templates,array,options)
if(mwa_dir == '/fred/oz048/MWA/'):
    generate_ozstar(infile,basename,n_subbands,rts_templates,array,options)

