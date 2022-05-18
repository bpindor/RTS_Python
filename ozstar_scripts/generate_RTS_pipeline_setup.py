def generate_galaxy(infile,basename,n_subbands,rts_templates,options):

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
        sourcelist_string =  '--sourcelist=/group/mwaeor/bpindor/PUMA/srclists/srclist_puma-v2_complete.txt'

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

    out_file = open('RTS_pipeline_setup_%d.sh' % options.chunk_number,'w+')

    out_file.write('#!/bin/bash -l\n')
    
    out_file.write('#SBATCH --job-name=RTS_setup\n')
    out_file.write('#SBATCH -o RTS_setup-%A.out\n')
    out_file.write('#SBATCH --ntasks=1\n')
    out_file.write('#SBATCH --ntasks-per-node=1\n')
    out_file.write('#SBATCH --time=00:60:00\n')
    out_file.write('#SBATCH --partition=gpuq\n')
    out_file.write('#SBATCH --account=mwaeor\n')
    out_file.write('#SBATCH --export=NONE\n')

    if(options.dynamic):
        out_file.write('aprun generate_dynamic_RTS_sourcelists.py  %s  -n %d --obslist=%s\n' % (sourcelist_string,options.dynamic,infile))

    # This is required to write auto_inner after data have been downloaded
    if(options.auto_gen):
        out_file.write('aprun generate_mwac_qRTS_auto.py %s %s %s %s --auto --chunk_number=%d %s %s %s %s %s %s %s %s %s\n' % (infile, basename, n_subbands, rts_templates, options.chunk_number,do_erase,flags_string,dynamic_string,tile_flags_string,dev_rts_string,tag_uvfits_string,tag_logs_string,sourcelist_string,process_string))

#    out_file.write('aprun generate_sReflag_mwaf.py %s \n' % infile)

#    out_file.write('mv sReflag_mwaf_files.sh sReflag_mwaf_files_%d.sh \n' % options.chunk_number)

#    out_file.write('aprun generate_sCompress_gpubox.py %s \n' % infile)

#    out_file.write('mv sCompress_gpubox_files.sh sCompress_gpubox_files_%d.sh\n' % options.chunk_number)

    out_file.close()

def generate_ozstar(infile,basename,n_subbands,rts_templates,options):

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
        sourcelist_string =  '--sourcelist=${MWA_DIR}/CODE/srclists/%s' % options.sourcelist
    else:
        sourcelist_string =  '--sourcelist=/group/mwaeor/bpindor/PUMA/srclists/srclist_puma-v2_complete.txt'

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

    n_obs = sum(1 for line in open(infile))
        
    out_file = open('RTS_pipeline_setup_%d.sh' % options.chunk_number,'w+')

    out_file.write('#!/bin/bash -l\n')
    
    out_file.write('#SBATCH --job-name=RTS_setup\n')
    out_file.write('#SBATCH -o RTS_setup-%A.out\n')
    out_file.write('#SBATCH --ntasks=1\n')
    out_file.write('#SBATCH --ntasks-per-node=1\n')
    out_file.write('#SBATCH --array=1-%s\n' % n_obs)
    out_file.write('#SBATCH --time=00:60:00\n')
    out_file.write('#SBATCH --partition=skylake\n')
    out_file.write('#SBATCH --account=oz048\n')
    out_file.write('#SBATCH --mem=6000\n')
    out_file.write('#SBATCH --export=NONE\n')
    out_file.write('#SBATCH --cpus-per-task=2\n')


    if(options.dynamic):
        out_file.write('cd $SLURM_SUBMIT_DIR\n')
        out_file.write('obsid_file="%s"\n' % infile)
        out_file.write('while read line\ndo\n   obs_array+=(${line})\ndone < ${obsid_file}\n')
#        obs_id = ${obs_array[${SLURM_ARRAY_TASK_ID}-1]}
#        metafile = '%sdata/%s/%s.metafits' % (mwa_dir,obs_id,obs_id)
#        if not (os.access(metafile,os.R_OK)):
#            metafile = '%sdata/%s/%s_metafits_ppds.fits' % (mwa_dir,obs_id,obs_id)
#            if not (os.access(metafile,os.R_OK)):
#                print 'ERROR: Cannot locate metafits file for %s' % obs_id
        out_file.write('cd %sdata/${obs_array[${SLURM_ARRAY_TASK_ID}-1]} \n' % (mwa_dir))
        metafile = '%sdata/${obs_array[${SLURM_ARRAY_TASK_ID}-1]}/${obs_array[${SLURM_ARRAY_TASK_ID}-1]}_metafits_ppds.fits' % mwa_dir
        # Check if target sourcelist exists
        srclist_head = (options.sourcelist).split('.')[0]
        
        out_file.write('srclist_file = %sdata/${obs_array[${SLURM_ARRAY_TASK_ID}-1]}/%s_${obs_array[${SLURM_ARRAY_TASK_ID}-1]}_patch%d.txt' % (mwa_dir,srclist_head,options.dynamic))
        out_file.write('if [ -f \"$srclist_file\" ]; then\n')
        out_file.write('    echo \"$srclist_file exists\"\n')
        out_file.write('else \n')
        out_file.write('    srun --export=ALL --mem=5000 python ${MWA_DIR}/CODE/srclists/srclist_by_beam.py -m %s -n %d -s ${MWA_DIR}/CODE/srclists/%s\n' % (metafile,options.dynamic,options.sourcelist))
        out_file.write('fi\n')
                       
        #out_file.write('srun --export=ALL --mem=5000 python ${MWA_DIR}/CODE/srclists/srclist_by_beam.py -m %s -n %d -s ${MWA_DIR}/CODE/srclists/%s\n' % (metafile,options.dynamic,options.sourcelist))
            
    out_file.close()

    

##########

import sys,os, glob
from optparse import OptionParser,OptionGroup

mwa_dir = os.getenv('MWA_DIR','/astro/mwaeor/MWA/')

usage = 'Usage: generate_RTS_pipeline_setup.py [text file of obsIDs] [basename] [# of subbands] [rts_templates]'

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

(options, args) = parser.parse_args()

infile = args[0]
basename = args[1]
n_subbands = args[2]
rts_templates = args[3]



try:
    in_file = open(infile)

except IOError, err:
    'Cannot open input file %s\n',str(infile)

mwa_dir = os.getenv('MWA_DIR','/astro/mwaeor/MWA/')

if(mwa_dir == '/astro/mwaeor/MWA/'):
    generate_galaxy(infile,basename,n_subbands,rts_templates,options)
if(mwa_dir == '/fred/oz048/MWA/'):
    generate_ozstar(infile,basename,n_subbands,rts_templates,options)
