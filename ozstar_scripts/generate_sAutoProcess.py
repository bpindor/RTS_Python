#! /usr/bin/env python

import sys,os, glob
from optparse import OptionParser,OptionGroup
from astropy.io import fits

def get_chips_band(infile):

    mwa_dir = os.getenv('MWA_DIR','/scratch/partner1019/MWA/')

    try:
        in_file = open(infile)

    except IOError, err:
        'Cannot open input file %s\n',str(infile)

    obsid = ((in_file.readline()).split())[0]

    in_file.close()

    meta_fits = fits.open(mwa_dir + 'data/' + obsid + '/' + obsid + '.metafits')

    freqcent = meta_fits[0].header['FREQCENT']

    if (freqcent < 170.0):
            band = 0 # Low Band
    else:
            band = 1 # High Band

    return band
    

###

def generate_galaxy(options):

    cwd = os.getcwd()

    # For now, we just specify a list of obsids

    if options.obslist is None:
        print 'At present, must specify list of obsids to be processed. Exiting...'
        sys.exit(1)

    n_obs = sum(1 for line in open('%s' % (options.obslist)))

    if options.download is True:
        download = 1
    else:
        download = 0

    if options.starttime is None:
        starttime = ''
    else:
        starttime = '-a %s' % options.starttime

    if options.do_erase is False:
        do_erase = ''
    else:
        do_erase = '--erase_GPUfiles'

    if options.chunksize is not None:
        chunksize = int(options.chunksize)
        num_lines = sum(1 for line in open('%s' % (options.obslist)))
        n_chunks = num_lines / chunksize
        if(num_lines % chunksize):
            n_chunks += 1
        for i in range(n_chunks):
            command = 'sed -n %d,%dp %s > obslist_temp_%d.dat' % (1+i*chunksize,(i+1)*(chunksize), options.obslist,i)
            os.system(command)
    else:
        n_chunks = 1
        command = 'cp  %s obslist_temp_%d.dat' % (options.obslist,0)
        os.system(command)

    if(options.channel_flags != mwa_dir + '/RTS/utils/flags/flagged_channels_default.txt'):
        flags_string = '--channel_flags=' + options.channel_flags
    else:
        flags_string = ''

    if(options.dynamic):
        dynamic_string = '--dynamic_sourcelist=%d' % options.dynamic
    else:
        dynamic_string = ''

    if options.tile_flags is not None:
        tile_flags_string = '--tile_flags=' + options.tile_flags
    else:
        tile_flags_string = ''

    if options.keep_qa_images:
        keep_qa_images_string = '--keep_qa_images'
    else:
        keep_qa_images_string = ''

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

    if(options.sourcelist):
        sourcelist_string = '--sourcelist=%s' % options.sourcelist
    else:
        sourcelist_string = ''

    if(options.process_time):
        process_string = '--process_time=%d' % options.process_time
    else:
        process_string = ''

    if(options.rts_only):
        rts_only_string = '--rts_only'
    else:
        rts_only_string = ''



    # Write qsub script

    cwd = os.getcwd()

    autofile = 'sAutoProcess_%s.sh' % options.label 

    auto_file = open(autofile,'w+')

    auto_file.write('#!/bin/bash\n')

    for i in range(n_chunks):
        auto_file.write('generate_getNGASData.py %s/obslist_temp_%d.dat %s %s --chunk_number=%d %s\n' % (cwd, i, options.tagname, options.templatefile,i,rts_only_string))
        # Write RTS processing scripts on galaxy
    #    if(options.dynamic):
    #        auto_file.write('generate_dynamic_RTS_sourcelists.py  %s  -n %d --obslist=%s/obslist_temp_%d.dat\n' % (sourcelist_string,options.dynamic,cwd,i))
        auto_file.write('generate_mwac_qRTS_auto.py %s/obslist_temp_%d.dat %s %d %s --auto --chunk_number=%d %s %s %s %s %s %s %s %s %s\n' % (cwd,i, options.n_bands,options.tagname, options.templatefile, i,do_erase,flags_string,dynamic_string,tile_flags_string,dev_rts_string,tag_uvfits_string,tag_logs_string,sourcelist_string,process_string))
        auto_file.write('generate_sReflag_mwaf.py %s/obslist_temp_%d.dat \n' % (cwd,i))

        auto_file.write('mv sReflag_mwaf_files.sh sReflag_mwaf_files_%d.sh \n' % i)

    #    auto_file.write('generate_sCompress_gpubox.py %s/obslist_temp_%d.dat \n' % (cwd,i))

    #    auto_file.write('mv sCompress_gpubox_files.sh sCompress_gpubox_files_%d.sh\n' % i)

        auto_file.write('generate_RTS_pipeline_setup.py %s/obslist_temp_%d.dat %s %d %s --auto --chunk_number=%d %s %s %s %s %s %s %s %s %s\n'  % (cwd, i, options.tagname, options.n_bands,options.templatefile,i,do_erase,flags_string,dynamic_string,tile_flags_string,dev_rts_string,tag_uvfits_string,tag_logs_string,sourcelist_string,process_string)) 
        # Write script which submits RTS processing scripts from within getNGASdata script on zeus
        if(i==n_chunks-1):
            auto_file.write('generate_autoProcess_pipeline_depend.py --obslist=%s/obslist_temp_%d.dat --chunk_number=%d %s --do_QA --full_obslist=%s --tagname=%s\n' % (cwd, i,i,rts_only_string,options.obslist,options.tagname)) 
        else:
            auto_file.write('generate_autoProcess_pipeline_depend.py --obslist=%s/obslist_temp_%d.dat --chunk_number=%d %s\n' % (cwd, i,i,rts_only_string)) 

        auto_file.write('chmod +x autoProcess_pipeline_%d.sh\n' % i)

        auto_file.write('generate_mwac_qRTS_auto.py %s/obslist_temp_%d.dat %s %d %s --auto --chunk_number=%d %s %s %s %s %s %s %s %s %s\n' % (cwd, i, options.tagname, options.n_bands,options.templatefile,i,do_erase,flags_string,dynamic_string,tile_flags_string,dev_rts_string,tag_uvfits_string,tag_logs_string,sourcelist_string,process_string))

        auto_file.write('chmod +x sRTS_auto_inner_%d.sh\n' % i)
        auto_file.write('sbatch sGetNGASData_pipeline_%d.sh\n' % i)

        # End of Chunks



    auto_file.close()
    print 'Wrote script file: sAutoProcess_%s.sh' % options.label  

    cmd = 'chmod +x sAutoProcess_%s.sh' % options.label
    os.system(cmd)

def generate_ozstar(options):

    cwd = os.getcwd()

    # For now, we just specify a list of obsids

    if options.obslist is None:
        print 'At present, must specify list of obsids to be processed. Exiting...'
        sys.exit(1)

    n_obs = sum(1 for line in open('%s' % (options.obslist)))

    if options.download is True:
        download = 1
    else:
        download = 0

    if options.starttime is None:
        starttime = ''
    else:
        starttime = '-a %s' % options.starttime

    if options.do_erase is False:
        do_erase = ''
    else:
        do_erase = '--erase_GPUfiles'

    if options.chunksize is not None:
        chunksize = int(options.chunksize)
        num_lines = sum(1 for line in open('%s' % (options.obslist)))
        n_chunks = num_lines / chunksize
        if(num_lines % chunksize):
            n_chunks += 1
        for i in range(n_chunks):
            command = 'sed -n %d,%dp %s > obslist_temp_%d.dat' % (1+i*chunksize,(i+1)*(chunksize), options.obslist,i)
            os.system(command)
    else:
        n_chunks = 1
        command = 'cp  %s obslist_temp_%d.dat' % (options.obslist,0)
        os.system(command)

    if(options.channel_flags != mwa_dir + '/RTS/utils/flags/flagged_channels_default.txt'):
        flags_string = '--channel_flags=' + options.channel_flags
    else:
        flags_string = ''

    if(options.dynamic):
        dynamic_string = '--dynamic_sourcelist=%d' % options.dynamic
    else:
        dynamic_string = ''

    if options.tile_flags is not None:
        tile_flags_string = '--tile_flags=' + options.tile_flags
    else:
        tile_flags_string = ''

    if options.keep_qa_images:
        keep_qa_images_string = '--keep_qa_images'
    else:
        keep_qa_images_string = ''

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

    if(options.sourcelist):
        sourcelist_string = '--sourcelist=%s' % options.sourcelist
    else:
        sourcelist_string = ''

    if(options.process_time):
        process_string = '--process_time=%d' % options.process_time
    else:
        process_string = ''

    if(options.rts_only):
        rts_only_string = '--rts_only'
    else:
        rts_only_string = ''   

    # run sbatch script
        
    autofile = 'sAutoProcess_%s.sh' % options.label 

    auto_file = open(autofile,'w+')

    auto_file.write('#!/bin/bash\n')

    for i in range(n_chunks):
        
        auto_file.write('generate_getGPUBOX_Data.py %s/obslist_temp_%d.dat --chunk_number=%d \n' % (cwd,i,i))
#        auto_file.write('generate_getNGASData.py %s/obslist_temp_%d.dat %s %s --chunk_number=%d %s\n' % (cwd, i, options.tagname, options.templatefile,i,rts_only_string))
        auto_file.write('generate_sRTS_auto.py %s/obslist_temp_%d.dat %s %d %s --auto --chunk_number=%d %s %s %s %s %s %s %s %s %s\n' % (cwd,i, options.tagname, options.n_bands,options.templatefile, i,do_erase,flags_string,dynamic_string,tile_flags_string,dev_rts_string,tag_uvfits_string,tag_logs_string,sourcelist_string,process_string))
#        auto_file.write('generate_mwac_qRTS_auto.py %s/obslist_temp_%d.dat %s %d %s --auto --chunk_number=%d %s %s %s %s %s %s %s %s %s\n' % (cwd,i, options.tagname, options.n_bands,options.templatefile, i,do_erase,flags_string,dynamic_string,tile_flags_string,dev_rts_string,tag_uvfits_string,tag_logs_string,sourcelist_string,process_string))
        auto_file.write('generate_sReflag_mwaf.py %s/obslist_temp_%d.dat \n' % (cwd,i))
        auto_file.write('mv sReflag_mwaf_files.sh sReflag_mwaf_files_%d.sh \n' % i)
        auto_file.write('generate_RTS_sourcelists.py %s/obslist_temp_%d.dat %s %d %s --auto --chunk_number=%d %s %s %s %s %s %s %s %s %s\n'  % (cwd, i, options.tagname, options.n_bands,options.templatefile,i,do_erase,flags_string,dynamic_string,tile_flags_string,dev_rts_string,tag_uvfits_string,tag_logs_string,sourcelist_string,process_string)) 
        #        auto_file.write('generate_RTS_pipeline_setup.py %s/obslist_temp_%d.dat %s %d %s --auto --chunk_number=%d %s %s %s %s %s %s %s %s %s\n'  % (cwd, i, options.tagname, options.n_bands,options.templatefile,i,do_erase,flags_string,dynamic_string,tile_flags_string,dev_rts_string,tag_uvfits_string,tag_logs_string,sourcelist_string,process_string)) 

#        if(i==n_chunks-1):
#            auto_file.write('generate_autoProcess_pipeline_depend.py --obslist=%s/obslist_temp_%d.dat --chunk_number=%d %s --do_QA --full_obslist=%s --tagname=%s --n_bands=%d\n' % (cwd, i,i,rts_only_string,options.obslist,options.tagname,options.n_bands)) 
#        else:
#            auto_file.write('generate_autoProcess_pipeline_depend.py --obslist=%s/obslist_temp_%d.dat --chunk_number=%d %s --n_bands=%d\n' % (cwd, i,i,rts_only_string,options.n_bands)) 

        auto_file.write('generate_autoProcess_pipeline_depend.py --obslist=%s/obslist_temp_%d.dat --chunk_number=%d %s --do_QA --full_obslist=%s --tagname=%s --n_bands=%d\n' % (cwd, i,i,rts_only_string,options.obslist,options.tagname,options.n_bands)) 
            
        auto_file.write('chmod +x autoProcess_pipeline_%d.sh\n' % i)

        auto_file.write('chmod +x sRTS_auto_inner_%d.sh\n' % i)
        auto_file.write('sbatch sgetGPUBOXData_pipeline_%d.sh\n' % i)

        # End of Chunks

    auto_file.close()
    print 'Wrote script file: sAutoProcess_%s.sh' % options.label  

    cmd = 'chmod +x sAutoProcess_%s.sh' % options.label
    os.system(cmd)
        


        
""" 
generate_sAutoProcess.py
Generates a sbatch script to automatically download and process a set of EOR observations.
"""

import sys,os, glob
from optparse import OptionParser,OptionGroup

mwa_dir = os.getenv('MWA_DIR','/scratch/partner1019/MWA/')

usage = 'Usage: generate_sAutoProcess.py [options]\n'
usage += '\tGenerates a sbatch script to automatically download\n'
usage += '\tand process a set of MWA observations.\n'

parser = OptionParser(usage=usage)
parser.add_option('-d','--label',dest='label',
                      help='Optional Label to be applied to generated script [default=%default]',type='string', default='Test')
parser.add_option('--template', dest='templatefile',
                      help='List of RTS templates',type='string')
parser.add_option('--tagname', dest='tagname',default='autoCals',
                      help='Tag string used to identify this processing run [default=%default]')
parser.add_option('--no_download',dest='download',default=True,
                      action='store_false',
                      help='Download GPU files from NGAS [default=%default]')
parser.add_option('--starttime',dest='starttime',default=None,
                      help='Time to start qsub job. The default is to start immediately')
parser.add_option('--erase_GPUfiles',action='store_true',dest='do_erase',default=False,help='Erase raw GPU files after processing [default=%default]')

parser.add_option('--obslist',dest='obslist',default=None,
                      help='List of EOR Obsids to Process')
parser.add_option('--chunksize',dest='chunksize',default=None,
                      help='Number of Obsids to process within each job')
parser.add_option('--disable_QA',dest='do_QA',action='store_false',default=True,help='Run Miriad Based QA on uvdump files')
parser.add_option('--run_CHIPS',dest='run_CHIPS',action='store_true',default=False,help='Run CHIPS PS Module')
parser.add_option('--channel_flags',dest="channel_flags",default=mwa_dir + '/RTS/utils/flags/flagged_channels_default.txt',
                      help="File containing flagged channels",metavar="CHANNELFLAGS") 
parser.add_option('--dynamic_sourcelist',dest='dynamic',type='int',default=0,help='Use Dynamically Generated Sourcelist for Calibration [number of sources] (only applies to first entry in template file)')                 
parser.add_option('--sourcelist',dest='sourcelist',type='string',default='',help='Specify the base catalog to be used to generate the dynamic sourcelist (only applies to first entry in template file) default is /astro/mwaeor/bpindor/PUMA/srclists/srclist_puma-v2_complete.txt')
parser.add_option('--tile_flags',dest="tile_flags",default=None,
                      help="Comma separated list of additional tiles to flag. Use tile names e.g. 8th tile on Rx10 is Tile108, so --tile_flags=108. Default is none.")
parser.add_option('--keep_qa_images',dest="keep_qa_images",default=False,action='store_true',
                      help="QA images are exported to FITS and not deleted. Default=%default.")
parser.add_option('--dev_rts',dest="dev_rts",default=False,action='store_true',
                  help="Use development version of the RTS")
parser.add_option('--tag_uvfits',dest="tag_uvfits",default=False,action='store_true',
                  help="Tag and copy uvfits to data subdirectory")
parser.add_option('--tag_logs',dest="tag_logs",default=False,action='store_true',
                  help="Tag and copy uvfits to data subdirectory")
parser.add_option('--reflag_mwaf',dest="reflag_mwaf",default=False,action='store_true',
                  help="Run reflagging to identify fine channels with large flagged fraction in Cotter mwaf files")
parser.add_option('--process_time',dest="process_time",type='int',default=0,help='Processing time per RTS run per obsid. Default is 20 minutes.')
parser.add_option('--band',dest='band',type='int',default=-1,help='EOR band for running CHIPS [0=low band, 1 = high band]')
parser.add_option('--nofilter_for_chips',dest='nofilter_for_chips',action='store_true',default=False,help='Run CHIPS over all the obsdis without qa_check')
parser.add_option('--email',dest='email_addr',default=None,help='Send notification email when jobs are completed/stopped')
parser.add_option('--rts_only',dest="rts_only",default=False,action='store_true',help="Skip downloading and RTS setup steps")
parser.add_option('--n_bands',dest='n_bands',type='int',default=24,help='Number of Coarse Bands')

(options, args) = parser.parse_args()

mwa_dir = os.getenv('MWA_DIR','/scratch/partner1019/MWA/')

if(mwa_dir == '/astro/mwaeor/MWA/'):
    generate_galaxy(options)
if(mwa_dir == '/fred/oz048/MWA/'):
    generate_ozstar(options)
