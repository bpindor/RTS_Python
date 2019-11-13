#! /usr/bin/env python

import sys,os,glob
import astropy.io.fits as pyfits
import subprocess
from optparse import OptionParser,OptionGroup

def match_cal_name(fieldname,sourceCat):

	rts_dir = os.getenv('RTSDIR','/scratch/astronomy556/MWA/RTS/')
	cals_file = open(rts_dir + 'utils/source_lists/MWA_Cals_Match.txt')

	objname = fieldname.split('_')[1]

	cal_name = None

	for line in cals_file:
		if(line.split()[0] == objname):
			if(sourceCat == 'srclist_culgoora_mod.txt'):
				cal_name = line.split()[3]
			if(sourceCat == 'pks.txt'):
				cal_name =  line.split()[1]

	return cal_name

#function that gives you the 'rts tile number' from the tile name
def tilename_to_rtstile(tilename):
   tiles = {
      'Tile104':0,
      'Tile103':1,
      'Tile102':2,
      'Tile101':3,  
      'Tile108':4,
      'Tile107':5,
      'Tile106':6,
      'Tile105':7,
      'Tile074':8,
      'Tile073':9,
      'Tile072':10,
      'Tile071':11,
      'Tile078':12,
      'Tile077':13,
      'Tile076':14,
      'Tile075':15,
      'Tile164':16,
      'Tile163':17,
      'Tile162':18,
      'Tile161':19,
      'Tile168':20,
      'Tile167':21,
      'Tile166':22,
      'Tile165':23,
      'Tile154':24,
      'Tile153':25,
      'Tile152':26,
      'Tile151':27,
      'Tile158':28,
      'Tile157':29,
      'Tile156':30,
      'Tile155':31,
      'Tile024':32,
      'Tile023':33,
      'Tile022':34,
      'Tile021':35,
      'Tile028':36,
      'Tile027':37,
      'Tile026':38,
      'Tile025':39,
      'Tile014':40,
      'Tile013':41,
      'Tile012':42,
      'Tile011':43,
      'Tile018':44,
      'Tile017':45,
      'Tile016':46,
      'Tile015':47,
      'Tile094':48,
      'Tile093':49,
      'Tile092':50,
      'Tile091':51,
      'Tile098':52,
      'Tile097':53,
      'Tile096':54,
      'Tile095':55,
      'Tile084':56,
      'Tile083':57,
      'Tile082':58,
      'Tile081':59,
      'Tile088':60,
      'Tile087':61,
      'Tile086':62,
      'Tile085':63,
      'Tile124':64,
      'Tile123':65,
      'Tile122':66,
      'Tile121':67,
      'Tile128':68,
      'Tile127':69,
      'Tile126':70,
      'Tile125':71,
      'Tile114':72,
      'Tile113':73,
      'Tile112':74,
      'Tile111':75,
      'Tile118':76,
      'Tile117':77,
      'Tile116':78,
      'Tile115':79,
      'Tile144':80,
      'Tile143':81,
      'Tile142':82,
      'Tile141':83,
      'Tile148':84,
      'Tile147':85,
      'Tile146':86,
      'Tile145':87,
      'Tile134':88,
      'Tile133':89,
      'Tile132':90,
      'Tile131':91,
      'Tile138':92,
      'Tile137':93,
      'Tile136':94,
      'Tile135':95,
      'Tile034':96,
      'Tile033':97,
      'Tile032':98,
      'Tile031':99,
      'Tile038':100,
      'Tile037':101,
      'Tile036':102,
      'Tile035':103,
      'Tile044':104,
      'Tile043':105,
      'Tile042':106,
      'Tile041':107,
      'Tile048':108,
      'Tile047':109,
      'Tile046':110,
      'Tile045':111,
      'Tile064':112,
      'Tile063':113,
      'Tile062':114,
      'Tile061':115,
      'Tile068':116,
      'Tile067':117,
      'Tile066':118,
      'Tile065':119,
      'Tile054':120,
      'Tile053':121,
      'Tile052':122,
      'Tile051':123,
      'Tile058':124,
      'Tile057':125,
      'Tile056':126,
      'Tile055':127
   }
   return tiles[tilename]

"""
generate_RTS_in.py
Takes as input header.txt file created by convert_ngas and generates rts .in files 
"""

import sys,os,glob
import subprocess
from optparse import OptionParser,OptionGroup
from math import sin, cos, atan2, asin, degrees, radians 

mwa_dir = os.getenv('MWA_DIR','/scratch/astronomy556/MWA/')

usage= 'Usage: generate_RTS_in_auto.py [data_dir] [input_file] [basename] [# of subbands] [array]\n'
usage+= 'Generate RTS .in files for use in Fornax/Galaxy RTS pipeline\n'

parser = OptionParser(usage=usage)

parser.add_option('--header',dest="input_file",default='header.txt',
                      help="corr2uvfits style header file",metavar="DATADIR")
parser.add_option('--templates',dest="template_list",default=None,
                      help="List of RTS template file",metavar="RTSLIST")
parser.add_option('--channel_flags',dest="channel_flags",default=mwa_dir + '/RTS/utils/flags/flagged_channels_default.txt',
                      help="File containing flagged channels",metavar="CHANNELFLAGS")
parser.add_option('--tile_flags',dest="tile_flags",default=None,
                      help="Comma separated list of additional tiles to flag. Use tile names e.g. 8th tile on Rx10 is Tile108, so --tile_flags=108. Default is none.")

(options, args) = parser.parse_args()


data_dir = args[0]
basename = args[1]
n_subbands = int(args[2])
array = str(args[3])

half_channelwidth = 0.02 # 20 khz

gps = data_dir[len(mwa_dir + 'data/'):]
# strip trailing slash if necessary
if(gps[-1] == '/'):
   gps = gps[:-1]

if not (array == '32T' or array == '128T'):
        print 'Array Parameter should be \'32T\' or \'128T\''
        exit(1)

# Read in list of RTS template files

if options.template_list is None:
    template_files = [mwa_dir + 'RTS/utils/templates/RTS_template_regrid.in']
else:
    try: 
        template_list_file = open(options.template_list)
    except IOError, err:
        'Cannot open list of RTS template files %s\n' % str(options.template_list)
        
    template_files = []
    for line in template_list_file:
        template_files.append(line.strip())
    template_list_file.close()

# Set channel flag file
if (options.channel_flags is not None):
   cmd = 'cp %s flagged_channels.txt' % (options.channel_flags)
   os.system(cmd)

#Generate tile flag file
try:
        in_file = pyfits.open(options.input_file)
except IOError, err:
        'Cannot open metadata file %s\n' % str(options.input_file)

rts_inputs = in_file[1].data['Input']
tile_names = in_file[1].data['Tile']
rts_flags = rts_inputs/2
tiles2rts = dict(zip(tile_names,rts_flags))



if (options.tile_flags is not None):
   cmd = 'rm flagged_tiles.txt'
   os.system(cmd)
   tile_flag_string=options.tile_flags
   flagged_tile_list=tile_flag_string.split(',')
   for tile in flagged_tile_list:
      tile_name='Tile'+'{0:03d}'.format(int(tile))
      #      rts_tile=str(tilename_to_rtstile(tile_name))
      rts_tile = str(tiles2rts[int(tile)])
      cmd = 'echo "%s" >> flagged_tiles.txt' % (rts_tile)
      os.system(cmd)

# Read in corr2uvfits style header file

if(options.input_file=='header.txt'):

	try:
		in_file = open(options.input_file)
	except IOError, err:
		'Cannot open header file %s\n' % str(options.input_file)

	for line in in_file:
	     if(line.find('N_CHANS')==0):
		 n_chans = int((line.split())[1])
	     if(line.find('N_SCANS')==0):
		 n_scans = int((line.split())[1])
	     if(line.find('INT_TIME')==0):
		 int_time = float((line.split())[1])
	     if(line.find('BANDWIDTH')==0):
		 full_bandwidth = float((line.split())[1])
		 bandwidth = full_bandwidth / float(n_subbands)
		 half_channelwidth = (full_bandwidth / (n_chans * 2.0)) / 1000.0 
	     if(line.find('FREQCENT')==0):
		 freq_cent = float((line.split())[1])
	     if(line.find('HA_HRS')==0):
		 ha_point = float((line.split())[1])
		 if(ha_point > 12.0):
		     ha_point -= 24.0
	     if(line.find('DEC_DEGS')==0):
		 dec_phase = float((line.split())[1])
		 dec_point = dec_phase
	     if(line.find('RA_HRS')==0):
		 ra_phase = float((line.split())[1])
	     if(line.find('FIELDNAME')==0):
		 fieldname = ((line.split())[1])


	in_file.close()
	
else:
	try:
		in_file = pyfits.open(options.input_file)
	except IOError, err:
		'Cannot open metadata file %s\n' % str(options.input_file)

	n_chans = int((in_file[0].header)['NCHANS'])
	n_scans = int((in_file[0].header)['NSCANS'])
	int_time = float((in_file[0].header)['INTTIME'])
	full_bandwidth = float((in_file[0].header)['BANDWDTH'])
	bandwidth = full_bandwidth / float(n_subbands)
	freq_cent = float((in_file[0].header)['FREQCENT'])
	ra_point = float((in_file[0].header)['RA'])
	dec_point = float((in_file[0].header)['DEC'])
        try:
            ra_phase=in_file[0].header['RAPHASE']
        except:
            print 'No RAPHASE'
            ra_phase = ra_point
            pass
	#ra_phase = float((in_file[0].header)['RAPHASE'])
	ra_phase = ra_phase / 15.0 # Convert to hrs
        try:
            dec_phase = float((in_file[0].header)['DECPHASE'])
        except:
            print 'No DECPHASE'
            dec_phase = dec_point
            pass        
	#dec_phase = float((in_file[0].header)['DECPHASE'])
	lst = float((in_file[0].header)['LST'])
	ha_point = (lst - ra_point) / 15.0
	fieldname = (in_file[0].header)['FILENAME']
	finechan = float((in_file[0].header)['FINECHAN'])
	half_channelwidth = (finechan / 2.0) / 1000.0
	az = radians(float((in_file[0].header)['AZIMUTH']))
	el = radians(float((in_file[0].header)['ALTITUDE'])) 
	lat = -26.70331940 # Check this??
	lat = radians(lat)

	ha_point = degrees(atan2( -sin(az)*cos(el), cos(lat)*sin(el)-sin(lat)*cos(el)*cos(az))) / 15.0
        dec_point = degrees(asin( sin(lat)*sin(el)+cos(lat)*cos(el)*cos(az) ))

	in_file.close()
	
	

##


rts_file_index = 0

# How many GPU boxes are there?

gpu_files = glob.glob(data_dir + '/*gpubox*00.fits')  

if(len(gpu_files) > 0):
	#band_list = [int(filename[filename.find('gpubox')+6:-8]) for filename in gpu_files] 
	band_list = [int(filename[-10:-8]) for filename in gpu_files]

else:
	uvfile_list = glob.glob(data_dir + '/*uvfits')
	band_list = [int(filename[-9:-7]) for filename in uvfile_list]   

band_list.sort()

subband_string = ''

for band in band_list:
    subband_string = subband_string + str(band) + ','

subband_string = subband_string[:-1]

pwrs2 = [2**n for n in range(10)]

for file in template_files:

    template_file = open(file)

    outfile = basename + ('_rts_%d.in' % rts_file_index)
    rts_file_index += 1

    out_file = open(outfile,"w+")
    set_subbands = 0

    # Check to make sure the template is of mwac type

    is_mwac = 0
    gps = data_dir[len(mwa_dir + 'data/'):]
    # strip trailing slash if necessary
    if(gps[-1] == '/'):
	    gps = gps[:-1]
    times = os.popen('timeconvert.py --gps=%s' % gps).read()
    jd = (times.split('\n')[5]).split()[1]
				    

    # Make first pass through template_file to:
    # i) work out imaging cadence
    # ii) check if we need to see the primary calibrator
    # iii) Are we reading directly from gpubox files
    # iv) Are we using Cotter flags

    CorrDumpTime = CorrDumpsPerCadence = -1.0
    usePrimaryCalibrator = False
    readGpuboxDirect = False
    ImportCotterFlags = False

    for line in template_file:
	if(line.find('CorrDumpTime')==0):
		CorrDumpTime = float(line[len('CorrDumpTime='):])
	if(line.find('CorrDumpsPerCadence')==0):    
		CorrDumpsPerCadence = float(line[len('CorrDumpsPerCadence='):])
	if(line.find('PrimaryCalibrator')==0):
		usePrimaryCalibrator=True
	if(line.find('SourceCatalogueFile')==0):
		SourceCat = line[line.rfind('/')+1:-1]
	if(line.find('ReadGpuboxDirect=1')==0):
		readGpuboxDirect = True
        if(line.find('ImportCotterFlags=1')==0):
                ImportCotterFlags = True
        # NumberofIterations is now the adjustable parameter and
        # all other cadence values are populated accordingly
        if(line.find('NumberOfIterations')==0):
		n_iterations = int(line[len('NumberOfIterations='):])

    if(CorrDumpTime > 0 and CorrDumpsPerCadence):
            imaging_cadence = int_time * CorrDumpsPerCadence
#	    imaging_cadence = CorrDumpTime * CorrDumpsPerCadence
    else:
	    imaging_cadence = 8.0 #default

    scan_time = n_scans * int_time

    if(n_iterations == 1 and (n_scans / n_iterations) > 64):
	    CorrDumpsPerCadence = 128
    else:
	    # Make sure number of CorrDumpsPerCadence is a power of 2
	    for p in pwrs2:
		    if(n_scans / n_iterations >= p):
			    CorrDumpsPerCadence = p
			    
	    #CorrDumpsPerCadence = n_scans / n_iterations

    
#    n_iterations = int(scan_time / imaging_cadence)

    ## HERE ##

#    print 'scan_time', scan_time
#    print 'n_iterations', n_iterations
#    print 'imaging_cadence', imaging_cadence

#    if(usePrimaryCalibrator):
#	    PrimaryCalibrator = match_cal_name(fieldname,SourceCat)
	    
    template_file.close()


#If using Cotter files, check if Cotter flag file exists
#Check for .mwaf files and only unzip if not already present.

    if ImportCotterFlags:
       mwaf_files_exist=False
       if (os.path.isfile(data_dir+'/'+gps+'_01.mwaf') and os.path.isfile(data_dir+'/'+gps+'_02.mwaf') and os.path.isfile(data_dir+'/'+gps+'_03.mwaf') and os.path.isfile(data_dir+'/'+gps+'_04.mwaf') and os.path.isfile(data_dir+'/'+gps+'_05.mwaf') and os.path.isfile(data_dir+'/'+gps+'_06.mwaf') and os.path.isfile(data_dir+'/'+gps+'_07.mwaf') and os.path.isfile(data_dir+'/'+gps+'_08.mwaf') and os.path.isfile(data_dir+'/'+gps+'_09.mwaf') and os.path.isfile(data_dir+'/'+gps+'_10.mwaf') and os.path.isfile(data_dir+'/'+gps+'_11.mwaf') and os.path.isfile(data_dir+'/'+gps+'_12.mwaf') and os.path.isfile(data_dir+'/'+gps+'_13.mwaf') and os.path.isfile(data_dir+'/'+gps+'_14.mwaf') and os.path.isfile(data_dir+'/'+gps+'_15.mwaf') and os.path.isfile(data_dir+'/'+gps+'_16.mwaf') and os.path.isfile(data_dir+'/'+gps+'_17.mwaf') and os.path.isfile(data_dir+'/'+gps+'_18.mwaf') and os.path.isfile(data_dir+'/'+gps+'_19.mwaf') and os.path.isfile(data_dir+'/'+gps+'_20.mwaf') and os.path.isfile(data_dir+'/'+gps+'_21.mwaf') and os.path.isfile(data_dir+'/'+gps+'_22.mwaf') and os.path.isfile(data_dir+'/'+gps+'_23.mwaf') and os.path.isfile(data_dir+'/'+gps+'_24.mwaf')):
          mwaf_files_exist=True
       if mwaf_files_exist is not True:
          cotter_flag_file_exists=False
          cotter_zip_file1=data_dir+'/'+gps+'_flags.zip'
          cotter_zip_file2=data_dir+'/'+gps+'/'+gps+'_flags.zip'
          if (os.path.exists(cotter_zip_file1)):
             cotter_flag_file_exists=True
             print 'unzipping Cotter flag file from %s.' % (cotter_zip_file1)
             cmd = 'unzip -o -j %s -d %s' % (cotter_zip_file1,data_dir)
             os.system(cmd)
             mwaf_files_exist=True
          elif (os.path.exists(cotter_zip_file2)):
             cotter_flag_file_exists=True
             print 'unzipping Cotter flag file from %s.' % (cotter_zip_file2)
             cmd = 'unzip -o -j %s -d %s' % (cotter_zip_file2,data_dir)
             os.system(cmd)
             mwaf_files_exist=True
          else:
             print 'Cotter flag zip file for obsid %s does not exist. WARNING: not using Cotter flags.' % (gps)



    # Now write new .in file

    template_file = open(file)

    for line in template_file:

        line_out = line
        if(line.find('BaseFilename')==0):
		line_out = 'BaseFilename=' + data_dir + '/*_gpubox' + '\n'
        if(line.find('NumberOfChannels')==0):
		line_out = line.replace(line[len('NumberOfChannels='):],str(n_chans/24) + '\n')
        if(line.find('ObservationFrequencyBase')==0):
            line_out = line.replace(line[len('ObservationFrequencyBase='):],str(freq_cent - full_bandwidth/2.0 - half_channelwidth) + '\n')
        if(line.find('ObservationPointCentreHA')==0):
            line_out = line.replace(line[len('ObservationPointCentreHA='):],str(ha_point) + '\n')
        if(line.find('ObservationPointCentreDec')==0):
            line_out = line.replace(line[len('ObservationPointCentreDec='):],str(dec_point) + '\n')
        if(line.find('ObservationImageCentreRA')==0):
            line_out = line.replace(line[len('ObservationImageCentreRA='):],str(ra_phase) + '\n')
        if(line.find('ObservationImageCentreDec')==0):
            line_out = line.replace(line[len('ObservationImageCentreDec='):],str(dec_phase) + '\n')
	
        if(line.find('NumberOfIterations')==0):
		if(array=='32T'):
			# CorrDumpsPerCadence is currently hard set to 4
			line_out = line.replace(line[len('NumberOfIterations='):],str(n_scans/4) + '\n')
		else:
			# Clip last iteration in case of missing data. TO DO: This is a hack.
			if(n_iterations > 2 and readGpuboxDirect==False):
				line_out = line.replace(line[len('NumberOfIterations='):],str(n_iterations-1) + '\n')
			else:
				line_out = line.replace(line[len('NumberOfIterations='):],str(n_iterations) + '\n')
        if(line.find('SubBandIDs=')==0):
		line_out = line.replace(line[len('SubBandIDs='):],subband_string + '\n')
		set_subbands = 1

	#MWAC

	if(line.find('UseCorrelatorInput=1')==0):
		is_mwac=1
	if(line.find('ObservationTimeBase=')==0):
		line_out = line.replace(line[len('ObservationTimeBase='):],jd + '\n')

		
	# Primary calibrator
        # This is now out of date and has been commented out

#	if(line.find('PrimaryCalibrator')==0):
#		line_out = line.replace(line[len('PrimaryCalibrator='):],PrimaryCalibrator + '\n')


	# Cadence

	if(line.find('CorrDumpTime')==0):
		line_out = line.replace(line[len('CorrDumpTime='):],str(int_time) + '\n')
        if(line.find('CorrDumpsPerCadence')==0):
		line_out = line.replace(line[len('CorrDumpsPerCadence='):],str(CorrDumpsPerCadence) + '\n')

        # Cotter Flags
        # If using cotter flags and flag files exist, turn off the RTS RFI flagging
        if(line.find('doRFIflagging=1')==0 and ImportCotterFlags and mwaf_files_exist):
                line_out = line.replace(line[len('doRFIflagging='):],str(0) + '\n')       
        if(line.find('ImportCotterFlags=1')==0 and ImportCotterFlags and mwaf_files_exist):
#                line_out = line.replace(line[len('ImportCotterFlags='):],str(1) + '\n'+'ImportCotterBasename='+data_dir+'/'+gps+'\n'+'//ReadMetafitsFile=1'+'\n'+'//MetafitsFilename='+data_dir+'/'+gps+'\n')
                line_out = line.replace(line[len('ImportCotterFlags='):],str(1) + '\n'+'ImportCotterBasename='+data_dir+'/RTS_'+gps+'\n'+'ReadMetafitsFile=1'+'\n'+'MetafitsFilename='+data_dir+'/'+gps+'\n')
        #if flag file does not exist set  ImportCotterFlags=0
        if(line.find('ImportCotterFlags=1')==0 and not mwaf_files_exist):
                #line_out = line.replace(line[len('ImportCotterFlags='):],str(0) +'\n')
                line_out = line.replace('ImportCotterFlags=1','//ImportCotterFlags=1\n')
	# Write out modified .in entry

        out_file.write(line_out)

    if(array=='128T'):
        if(set_subbands == 0):
            out_file.write('SubBandIDs='+subband_string+'\n')
    out_file.close()
    template_file.close()


            
            
            
