#! /usr/bin/env python

def generate_fornax(infile,rts_job):

    n_obs = sum(1 for line in open(infile))

    job_number = (rts_job.split('.'))[0]

    data_dir = mwa_dir + 'data/%s' % obs_id 
 
    out_file = open('qQA_pipeline.sh','w+')

    out_file.write('#!/bin/bash\n')
    out_file.write('#PBS -l select=1:ncpus=12:mem=8gb\n')
    out_file.write('#PBS -l walltime=00:%d:00\n' % (n_obs * 10))
    out_file.write('#PBS -m e\n')
    out_file.write('#PBS -q workq\n')
#    out_file.write('#PBS -q R2406045\n')
    out_file.write('#PBS -W group_list=partner1019\n')
    out_file.write('source /scratch2/mwaeor/bpindor/MWA_Python/bin/activate \n')
    out_file.write('cd $PBS_O_WORKDIR\n')

    out_file.write('/group/mwaops/CODE/bin/RTS_cal_plot.py --path=$MWA_DIR/data/ --list=%s --CheckFR --niter=7 --NumCal=300 --JobOut=$PBS_O_WORKDIR/qRTS_auto_wrapp.o%s\n' % (infile,job_number))
    out_file.write('pipeline_qa.py %s\n' % infile)
    out_file.write('/group/mwaops/CODE/bin/RTS_cal_plot.py --path=$MWA_DIR/data/ --list=%s --PlotCal\n' % infile)

    out_file.close()

def generate_galaxy(infile,rts_job,options):

    n_obs = sum(1 for line in open(infile))

    job_number = (rts_job.split('.'))[0]

    #batch file for imaging in parallel
    out_file = open('qQA_pipeline_image.sh','w+')

    out_file.write('#!/bin/bash -l\n')
    out_file.write('#SBATCH -J qa_image\n') 
    out_file.write('#SBATCH -o qa_image-%A-%a.out\n')
    out_file.write('#SBATCH --ntasks=1\n')
    out_file.write('#SBATCH --ntasks-per-node=1\n')
    out_file.write('#SBATCH --time=00:%d:00\n' % (30))
    out_file.write('#SBATCH --partition=gpuq\n')
    out_file.write('#SBATCH --account=mwaeor\n')
    out_file.write('#SBATCH --export=NONE\n')

#    out_file.write('source /scratch2/mwaeor/bpindor/MWA_Python/bin/activate \n')
    out_file.write('cd $SLURM_SUBMIT_DIR\n')

    out_file.write('obsid_file="%s"\n' % infile)
    out_file.write('while read line\ndo\n   obs_array+=(${line})\ndone < ${obsid_file}\n')
    out_file.write('pipeline_qa_uvfits.py ${obs_array[${SLURM_ARRAY_TASK_ID}-1]} '+ keep_qa_images_string+'\n')
#    out_file.write('pipeline_qa_images.py ${obs_array[${SLURM_ARRAY_TASK_ID}-1]} '+ keep_qa_images_string+'\n')


    out_file.close()

    #batch file for phase and amplitude plotting in parallel
    out_file = open('qQA_pipeline_phase_amp_plots.sh','w+')

    out_file.write('#!/bin/bash -l\n')
    out_file.write('#SBATCH -J qa_p_a_plots\n')
    out_file.write('#SBATCH -o qa_p_a_plots-%A-%a.out\n')
    out_file.write('#SBATCH --ntasks=1\n')
    out_file.write('#SBATCH --ntasks-per-node=1\n')
    out_file.write('#SBATCH --time=00:%d:00\n' % (10))
    out_file.write('#SBATCH --partition=gpuq\n')
    out_file.write('#SBATCH --account=mwaeor\n')
    out_file.write('#SBATCH --export=NONE\n')

#    out_file.write('source /scratch2/mwaeor/bpindor/MWA_Python/bin/activate \n')
    out_file.write('cd $SLURM_SUBMIT_DIR\n')

    out_file.write('obsid_file="%s"\n' % infile)
    out_file.write('while read line\ndo\n   obs_array+=(${line})\ndone < ${obsid_file}\n')
    #out_file.write('sql_make_metafiles.py -l --rts --gps=${obs_array[${SLURM_ARRAY_TASK_ID}-1]}  --header=${MWA_DIR}data/${obs_array[${SLURM_ARRAY_TASK_ID}-1]}/header.txt --antenna=antenna_locations.txt ${MWA_DIR}data/${obs_array[${SLURM_ARRAY_TASK_ID}-1]}/antenna_locations.txt --instr=${MWA_DIR}data/${obs_array[${SLURM_ARRAY_TASK_ID}-1]}/instr_config.txt\n')
    #out_file.write('aprun python /group/mwaeor/CODE/MWA_Tools/scripts/plot_CalSols.py  --path=$MWA_DIR/data --obsid=${obs_array[${SLURM_ARRAY_TASK_ID}-1]} --CheckFR --niter=7 --NumCal=300 --JobOut=$SLURM_SUBMIT_DIR/qRTS_auto_wrapp.o0\n')
    #out_file.write('aprun python /group/mwaeor/CODE/MWA_Tools/scripts/plot_CalSols.py --path=$MWA_DIR/data --obsid=${obs_array[${SLURM_ARRAY_TASK_ID}-1]} --PlotCal\n')
    out_file.write('aprun python /group/mwaeor/CODE/RTS/scripts/plot_CalSols.py  --base_dir=$MWA_DIR/data/${obs_array[${SLURM_ARRAY_TASK_ID}-1]} \n')

    #batch file for plotting rms and making pdf
    out_file = open('qQA_pipeline_plot_rms.sh','w+')

    out_file.write('#!/bin/bash -l\n')
    out_file.write('#SBATCH --ntasks=1\n')
    out_file.write('#SBATCH -o plot_rms-%A.out\n')
    out_file.write('#SBATCH --ntasks-per-node=1\n')
    out_file.write('#SBATCH --time=00:%d:00\n' % (n_obs * 1))
    out_file.write('#SBATCH --partition=gpuq\n')
    out_file.write('#SBATCH --account=mwaeor\n')
    out_file.write('#SBATCH --export=NONE\n')

#    out_file.write('source /scratch2/mwaeor/bpindor/MWA_Python/bin/activate \n')
    out_file.write('cd $SLURM_SUBMIT_DIR\n')

    out_file.write('obsid_file="%s"\n' % infile)
    out_file.write('aprun -n 1 pipeline_qa_plot_rms.py ${obsid_file}\n')
    out_file.write('cp qa_plots_pdf.pdf qa_plots_%s.pdf\n' % options.tagname)

    out_file.close()

def generate_ozstar(infile,options):
    
    n_obs = sum(1 for line in open(infile))

    job_number = (rts_job.split('.'))[0]

    #batch file for imaging in parallel
    out_file = open('qQA_pipeline_image.sh','w+')

    out_file.write('#!/bin/bash -l\n')
    out_file.write('#SBATCH -J qa_image\n') 
    out_file.write('#SBATCH -o qa_image-%A-%a.out\n')
    out_file.write('#SBATCH --ntasks=1\n')
    out_file.write('#SBATCH --ntasks-per-node=1\n')
    out_file.write('#SBATCH --time=00:%d:00\n' % (30))
    out_file.write('#SBATCH --partition=skylake\n')
    out_file.write('#SBATCH --account=oz048\n')
    out_file.write('#SBATCH --export=NONE\n')

    out_file.write('cd $SLURM_SUBMIT_DIR\n')

    out_file.write('obsid_file="%s"\n' % infile)
    out_file.write('while read line\ndo\n   obs_array+=(${line})\ndone < ${obsid_file}\n')
    out_file.write('srun --export=ALL /fred/oz048/bpindor/utils/pipeline_qa_uvfits.py ${obs_array[${SLURM_ARRAY_TASK_ID}-1]} --subdir=%s\n' % options.tagname)

    out_file.close()

    #batch file for phase and amplitude plotting in parallel
    out_file = open('qQA_pipeline_phase_amp_plots.sh','w+')

    out_file.write('#!/bin/bash -l\n')
    out_file.write('#SBATCH -J qa_p_a_plots\n')
    out_file.write('#SBATCH -o qa_p_a_plots-%A-%a.out\n')
    out_file.write('#SBATCH --ntasks=1\n')
    out_file.write('#SBATCH --ntasks-per-node=1\n')
    out_file.write('#SBATCH --time=00:%d:00\n' % (10))
    out_file.write('#SBATCH --partition=skylake\n')
    out_file.write('#SBATCH --account=oz048\n')
    out_file.write('#SBATCH --export=NONE\n')
    out_file.write('#SBATCH --mem=5000\n')

    out_file.write('cd $SLURM_SUBMIT_DIR\n')

    out_file.write('obsid_file="%s"\n' % infile)
    out_file.write('while read line\ndo\n   obs_array+=(${line})\ndone < ${obsid_file}\n')
    out_file.write('srun --export=ALL --mem=5000 python ${MWA_DIR}/CODE/RTS/scripts/plot_CalSols.py  --base_dir=$MWA_DIR/data/${obs_array[${SLURM_ARRAY_TASK_ID}-1]} --subdir=%s/ -i \n' % options.tagname)
    


    
#####

import sys,os, glob
from optparse import OptionParser,OptionGroup

usage = 'generate_qQA_pipeline.py [text file of obsIDs]'

parser = OptionParser(usage=usage)

parser.add_option('--keep_qa_images',action="store_true",dest="keep_qa_images", help="Export the qa images to FITS and keep them. Default is False (images are deleted)")
parser.add_option('--tagname', dest='tagname',default='autoCals',
                      help='Tag string used to identify this processing run [default=%default]')

(options, args) = parser.parse_args()

infile = args[0]
rts_job = args[1]

try:
    in_file = open(infile)

except IOError, err:
    'Cannot open input file %s\n',str(infile)

mwa_dir = os.getenv('MWA_DIR','/scratch/partner1019/MWA/')

if (options.keep_qa_images):
   keep_qa_images_string='--keep_qa_images'
else: 
   keep_qa_images_string=' '

if(mwa_dir == '/scratch/partner1019/MWA/'):
    generate_fornax(infile,rts_job,options)
if(mwa_dir == '/scratch2/mwaeor/MWA/'):
    generate_galaxy(infile,rts_job,options)
if(mwa_dir == '/fred/oz048/MWA/'):
    generate_ozstar(infile,options)


