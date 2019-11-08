#! /usr/bin/env python

def generate_galaxy(options):

    mwa_dir = os.getenv('MWA_DIR','/astro/mwaeor/MWA/')
    cwd = os.getcwd()

    out_file = open('autoProcess_pipeline_%d.sh' % options.chunk_number, 'w+')

    if(options.rts_only):
        out_file.write('RTS_JOB_%d=$(sbatch -M galaxy sRTS_auto_wrapper_%d.sh | cut -d \" \" -f 4)\n' % (options.chunk_number,options.chunk_number))
        out_file.write('echo RTS_JOB_%d:$RTS_JOB_%d\n' % (options.chunk_number,options.chunk_number))

    else:
        out_file.write('SETUP_JOB_%d=$(sbatch -M galaxy RTS_pipeline_setup_%d.sh | cut -d \" \" -f 4)\n' % (options.chunk_number,options.chunk_number))
        out_file.write('echo SETUP_%d:$SETUP_JOB_%d\n' % (options.chunk_number,options.chunk_number))
        out_file.write('REFLAG_JOB_%d=$(sbatch -M galaxy --dependency=afterok:$SETUP_JOB_%d sReflag_mwaf_files_%d.sh | cut -d \" \" -f 4)\n' % (options.chunk_number,options.chunk_number,options.chunk_number))
        out_file.write('echo REFLAG_JOB_%d:$REFLAG_JOB_%d\n' % (options.chunk_number,options.chunk_number))
        out_file.write('RTS_JOB_%d=$(sbatch -M galaxy --dependency=afterok:$REFLAG_JOB_%d sRTS_auto_wrapper_%d.sh | cut -d \" \" -f 4)\n' % (options.chunk_number,options.chunk_number,options.chunk_number))
        out_file.write('echo RTS_JOB_%d:$RTS_JOB_%d\n' % (options.chunk_number,options.chunk_number))

    #QA
    if(options.doQA):
        n_obs = sum(1 for line in open('%s' % (options.full_obslist)))
        out_file.write('generate_qQA_pipeline.py %s $RTS_JOB_%d --tagname=%s\n' % (options.full_obslist,options.chunk_number,options.tagname))
        out_file.write('QA_JOB1=$(sbatch -M galaxy --array=1-%s --dependency=afterok:$RTS_JOB_%d qQA_pipeline_image.sh | cut -d " " -f 4)\n' % (n_obs,options.chunk_number))
        out_file.write('echo $QA_JOB1\n')
        out_file.write('QA_JOB3=$(sbatch -M galaxy --array=1-%s --dependency=afterok:$RTS_JOB_%d qQA_pipeline_phase_amp_plots.sh | cut -d " " -f 4)\n' % (n_obs,options.chunk_number))
        out_file.write('echo $QA_JOB3\n')
        out_file.write('QA_JOB2=$(sbatch -M galaxy --dependency=afterany:$QA_JOB1 qQA_pipeline_plot_rms.sh | cut -d " " -f 4)\n')
        out_file.write('echo $QA_JOB2\n')
        out_file.write('LAST_JOB=$QA_JOB2\n')   

    out_file.close()

def generate_ozstar(options):

    mwa_dir = os.getenv('MWA_DIR','/astro/mwaeor/MWA/')
    cwd = os.getcwd()

    out_file = open('autoProcess_pipeline_%d.sh' % options.chunk_number, 'w+')

    if(options.rts_only):
        out_file.write('RTS_JOB_%d=$(sbatch  sRTS_auto_wrapper_%d.sh | cut -d \" \" -f 4)\n' % (options.chunk_number,options.chunk_number))
        out_file.write('echo RTS_JOB_%d:$RTS_JOB_%d\n' % (options.chunk_number,options.chunk_number))

    else:
        out_file.write('SETUP_JOB_%d=$(sbatch  RTS_pipeline_setup_%d.sh | cut -d \" \" -f 4)\n' % (options.chunk_number,options.chunk_number))
        out_file.write('echo SETUP_%d:$SETUP_JOB_%d\n' % (options.chunk_number,options.chunk_number))
        out_file.write('REFLAG_JOB_%d=$(sbatch  --dependency=afterok:$SETUP_JOB_%d sReflag_mwaf_files_%d.sh | cut -d \" \" -f 4)\n' % (options.chunk_number,options.chunk_number,options.chunk_number))
        out_file.write('echo REFLAG_JOB_%d:$REFLAG_JOB_%d\n' % (options.chunk_number,options.chunk_number))
        out_file.write('RTS_JOB_%d=$(sbatch  --dependency=afterok:$REFLAG_JOB_%d sRTS_auto_wrapper_%d.sh | cut -d \" \" -f 4)\n' % (options.chunk_number,options.chunk_number,options.chunk_number))
        out_file.write('echo RTS_JOB_%d:$RTS_JOB_%d\n' % (options.chunk_number,options.chunk_number))

    #QA
    if(options.doQA):
        
        #n_obs = sum(1 for line in open('%s' % (options.full_obslist)))
        #out_file.write('generate_sQA_pipeline.py %s $RTS_JOB_%d --tagname=%s --n_bands=%d\n' % (options.full_obslist,options.chunk_number,options.tagname,options.n_bands))
        # lets have each chunk do its own QA in case of failures or
        # in case where last chunk is not last to finish
        n_obs = sum(1 for line in open('%s' % (options.obslist)))
        out_file.write('generate_sQA_pipeline.py %s $RTS_JOB_%d --tagname=%s --n_bands=%d --chunk_number=%d\n' % (options.obslist,options.chunk_number,options.tagname,options.n_bands,options.chunk_number))
#        out_file.write('QA_JOB1=$(sbatch  --array=1-%s --dependency=afterok:$RTS_JOB_%d qQA_pipeline_image.sh | cut -d " " -f 4)\n' % (n_obs,options.chunk_number))
#        out_file.write('echo $QA_JOB1\n')
        out_file.write('QA_JOB1=$(sbatch  --array=1-%s --dependency=afterok:$RTS_JOB_%d sQA_pipeline_phase_amp_plots_%d.sh | cut -d " " -f 4)\n' % (n_obs,options.chunk_number,options.chunk_number))
        out_file.write('echo $QA_JOB1\n')
        out_file.write('QA_JOB2=$(sbatch  --array=1-%s --dependency=afterany:$RTS_JOB_%d sQA_pipeline_uvfits_rms_%d.sh | cut -d " " -f 4)\n' % (n_obs,options.chunk_number,options.chunk_number))
        out_file.write('echo $QA_JOB2\n')
        out_file.write('LAST_JOB=$QA_JOB2\n')   

    out_file.close()
    

##########

import sys,os, glob
from optparse import OptionParser,OptionGroup

mwa_dir = os.getenv('MWA_DIR','/astro/mwaeor/MWA/')

usage = 'Usage: generate_autoProcess_pipeline_depend.py [options]'

parser = OptionParser(usage=usage)

parser.add_option('--obslist',dest='obslist',default=None,
                      help='List of EOR Obsids to Process')
parser.add_option('--chunk_number',dest='chunk_number',type='int',default=0,
                  help='Chunk number for automatic processing')
parser.add_option('--tagname', dest='tagname',default='autoCals',
                      help='Tag string used to identify this processing run [default=%default]')
parser.add_option('--do_QA', action="store_true",dest='doQA',help="Run QA Jobs at the end of this (final) chunk")
parser.add_option('--full_obslist',dest='full_obslist',default=None,
                      help='Full list of EOR Obsids to Process (for QA)')
parser.add_option('--rts_only',dest="rts_only",default=False,action='store_true',help="Skip downloading and RTS setup steps")
parser.add_option('--n_bands',dest='n_bands',type='int',default=24,help='Number of Coarse Bands')


(options, args) = parser.parse_args()

infile = options.obslist

try:
    in_file = open(infile)

except IOError, err:
    'Cannot open input file %s\n',str(infile)

mwa_dir = os.getenv('MWA_DIR','/scratch/partner1019/MWA/')

if(mwa_dir == '/astro/mwaeor/MWA/'):
    generate_galaxy(options)
if(mwa_dir == '/fred/oz048/MWA/'):
    generate_ozstar(options)



    
    
    



    
        
        
