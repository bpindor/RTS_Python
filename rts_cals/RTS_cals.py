import numpy as np
from astropy.io import fits
import os

class rts_cal():
    def __init__(self,n_bands=24):
        self.n_ants = 128
        self.n_bands = n_bands
        self.antennas=[]
        self.flagged_antennas = []
        for i in range(self.n_ants):
            self.antennas.append(rts_antenna(i,n_bands=n_bands))
        self.freqs = None
        self.metafile = None
        self.BP_weights = [None] * n_bands

        

    # Each BP file contains all antennas, so load must occur at top level

    def load_BP_jones(self,band,path='./',raw=False,add_flagged=True,chan_bw=0.04):
        calsfile = path + '/BandpassCalibration_node%03d.dat' % band

        try :
            cals_file = open(calsfile)
        except IOError as e:
            raise

        present_freqs = ((cals_file.readline()).split(','))
        present_freqs = [int(round(float(p)/chan_bw)) for p in present_freqs]
        flagged_channels = np.array([i for i in range(32) if i not in present_freqs])
        all_gains = []

        all_antennas = []
        for line in cals_file:
            gains = line.split(',')
            antenna = int(gains[0])
            all_antennas.append(antenna-1)
            f_gains = [float(g) for g in gains]
            #amp_gains = np.array(map(float,gains[1::2]))
            #phase_gains = np.array(map(float,gains[2::2]))
            amp_gains = np.array(f_gains[1::2])
            phase_gains = np.array(f_gains[2::2])
            # These gains are in (surely) amp-phase, so let's convert to real,imag
            r_gains = amp_gains * np.cos(phase_gains)
            i_gains = amp_gains * np.sin(phase_gains)
            c_gains = r_gains + 1.0j * i_gains
            all_gains.append(c_gains)

        present_antennas = set(all_antennas)
        flagged_antennas = [a for a in range(self.n_ants) if a not in present_antennas]
        for f in flagged_antennas:
            if(f not in self.flagged_antennas):
                (self.flagged_antennas).append(f)

        # gains are in order (data, fits) for XX, XY, YX, YY

        if(raw==True):
            xx_gains = all_gains[::8]
            xy_gains = all_gains[2::8]
            yx_gains = all_gains[4::8]
            yy_gains = all_gains[6::8]
        else:
            xx_gains = [all_gains[i] for i in range(len(all_gains)) if i%8 == 1]
            xy_gains = [all_gains[i] for i in range(len(all_gains)) if i%8 == 3]
            yx_gains = [all_gains[i] for i in range(len(all_gains)) if i%8 == 5]
            yy_gains = [all_gains[i] for i in range(len(all_gains)) if i%8 == 7]

        antenna_number = [all_antennas[i] for i in range(len(all_antennas)) if i%8 == 0]

        # Need to treat separate frequency channels

        all_jones = []
        flagged_jones = np.array([[0,0],[0,0]])

        for i in range(len(xx_gains)):
            jones = [np.array([[xx,xy],[yx,yy]]) for xx,xy,yx,yy in zip(xx_gains[i],xy_gains[i],yx_gains[i],yy_gains[i])]
            self.antennas[antenna_number[i]].BP_jones[band-1] = jones
            if(add_flagged):
                # Add in flagged channels
                for ch in range(32):
                    if(ch in flagged_channels):
                        self.antennas[antenna_number[i]].BP_jones[band-1].insert(ch,flagged_jones)    
                self.BP_weights[band-1] = np.ones(32)
                self.BP_weights[band-1][flagged_channels] = 0.0
            else:
                self.BP_weights[band-1] = np.ones(len(present_freqs))


        cals_file.close()

    def load_DI_jones(self,band,path='./',filename='/DI_JonesMatrices_node'):
        calsfile = path + '%s%03d.dat' % (filename,band)

        try:
            cals_file = open(calsfile)
        except IOError as e:
            raise
            
        flux_scale = float(cals_file.readline())
        all_jones = []

        for line in cals_file:
            xxr, xxi, xyr, xyi, yxr, yxi, yyr, yyi = line.split(',')
            jones = np.array([[float(xxr) + 1.0j * float(xxi),float(xyr) + 1.0j * float(xyi)],[float(yxr) + 1.0j * float(yxi),float(yyr) + 1.0j * float(yyi)]]) 
            all_jones.append(jones)

        cals_file.close()

        DI_jones = [j for j in all_jones[1:]]
        DI_jones_ref = all_jones[0]

        for a,dj in zip(self.antennas,DI_jones):
            a.DI_jones[band-1] = dj
            a.DI_jones_ref[band-1] = DI_jones_ref

    def load_cals_from_single_jones(self,single_jones,antenna=0,cent_chan=12,mode='single_cal'):
        flagged_channels = [0,1,16,30,31]
        channels = np.arange(32)
        keep_channels = np.delete(channels,flagged_channels)
        ant = self.antennas[antenna]
        total_ch = 0
        base_ch = 0
        n_chan = 32
        band_0 = cent_chan - 12
        
        # single jones are in physical frequency order, so we have to 
        # go back to coarse PFB order for consistency with RTS cals
        for b in range(self.n_bands):
            band_number = b + cent_chan - 12
            if(band_number > 128):
                if(band_0 > 128):
                    band_index = self.n_bands -(b+1)
                else:
                    band_index = self.n_bands-(band_number-128)
            else:
                band_index = b
            if(mode == 'single_cal'):    
                ant.DI_jones_ref[band_index] = np.array([[1.0, 0.0],[0.0, 1.0]])
                ant.DI_jones[band_index] = np.array([[1.0, 0.0],[0.0, 1.0]])
                ant.BP_jones[band_index] = [None] * n_chan
                for ch_n,ch in enumerate(ant.BP_jones[band_index]):
                    ch_cals = [[[],[]],[[],[]]]
                    for pol1 in [0,1]:
                        for pol2 in [0,1]:
                            ch_cals[pol1][pol2] = single_jones[pol1][pol2][total_ch]
                            ant.BP_jones[band_index][ch_n] = np.array(ch_cals)
                    total_ch += 1
            else:

                # This should have already been set, so check DI_jones_ref exists
                # The calibration applied to data is
                # DI_Jones_ref * (DI_jones * BP_jones).I * V_ij
                # So (DI_ref).I * single_jones = (DI_Jones * BP_jones).I
                # single_jones.I * DI_ref = DI_jones * BP_jones
                # Average over frequency
                # <single_jones.I * DI_ref> ~ DI_jones
                # BP_jones = (DI_jones).I * (single_jones.I) * DI_ref

                if(ant.DI_jones_ref[band_index] is None):
                    print('Error: DI_jones_ref should be set but is None')

                #i) select out channels for this band
                
                single_jones_band = single_jones[base_ch:base_ch+n_chan,:,:]

                #ii) form single_jones.I * DI_ref

                avg_jones = np.zeros_like(single_jones_band)
                avg_jones[keep_channels] = np.matmul(np.linalg.inv(single_jones_band[keep_channels]),ant.DI_jones_ref[band_index])

                #iii) Average to get DI_jones

                ant.DI_jones[band_index] = np.mean(avg_jones[keep_channels,:,:],axis=0)
                ant.BP_jones[band_index] = np.zeros_like(single_jones_band)
                

                for ch in keep_channels:
                    ch_cals = np.matmul(np.linalg.inv(ant.DI_jones[band_index]),np.matmul(np.linalg.inv(single_jones[ch+base_ch,:,:]),ant.DI_jones_ref[band_index]))
                    ant.BP_jones[band_index][ch] = ch_cals
                
            base_ch += n_chan

    def average_BP_jones(self):
        for idx,a in enumerate(self.antennas):
            if(a!=None):
                self.antennas[idx].average_BP_jones()
        freqs = [(self.freqs[2*i] + self.freqs[2*i+1])/2.0 for i in range(len(self.freqs)/2)]
        self.freqs = freqs

    def form_single_jones(self,invert_order=False):
        for idx,a in enumerate(self.antennas):
            if(a!=None):
                self.antennas[idx].form_single_jones(invert_order=invert_order)

    def load_all_BP_jones(self,path='./',raw=False,add_flagged=True,n_bands=24):
        for i in range(1,n_bands+1):
            self.load_BP_jones(i,path=path,raw=raw,add_flagged=add_flagged)
            
    def load_all_DI_jones(self,path='./',filename='/DI_JonesMatrices_node',n_bands=24):
        for i in range(1,n_bands+1):
            self.load_DI_jones(i,path=path,filename=filename)

    def load_metadata(self,metafile):
        from astropy.time import Time
        self.metafile = metafile
        meta_file = fits.open(metafile)
        self.meta_gains = (meta_file[1].data)['Gains']
        t = Time(meta_file[0].header['DATE-OBS'])
        self.mjd = 2.4e6 + t.mjd
        # set frequencies across the whole band
        nchans = meta_file[0].header['nchans']
        fine_chan = meta_file[0].header['FINECHAN']
        freqcent = meta_file[0].header['FREQCENT']

        freqs = (fine_chan * 1e3) * np.arange(nchans)
        self.freqs = freqs + freqcent * 1e6 - (freqs[int(nchans/2)] + freqs[int(nchans/2) -1]) / 2.0 - (fine_chan * 1e3) / 2.0
        self.cent_chan = meta_file[0].header['CENTCHAN']
        meta_file.close()
 

    def apply_metafits_gains(self):
        for i in range(self.n_ants):
            g0 = float(self.meta_gains[2*i,0])
            norm_gains = [float(g)/g0 for g in self.meta_gains[2*i]]
            d = [di / g for di,g in zip(self.antennas[i].DI_jones,reversed(norm_gains))]
            self.antennas[i].DI_jones = d
        


    def all_BP_jones(self,antenna=0,reverse_bands=False,cent_chan=12,pol='xx'):
        if(pol=='xx'): pol_index=(0,0)
        if(pol=='xy'): pol_index=(0,1)
        if(pol=='yx'): pol_index=(1,0)
        if(pol=='yy'): pol_index=(1,1)
        
        all_BP_xx = []
        all_BP_yy = []
        a = self.antennas[antenna]
        band_0 = cent_chan - 12
        if(a!=None):
            for b,B in enumerate(a.BP_jones):
                band_number = b + cent_chan - 12
                if(band_number < 129):
                   if(B is not None):
                        for chan in B:
                            all_BP_xx.append(chan[0,0]) 
                else:
                    if(band_0 > 128):
                        #entire band is reversed
                        if(a.BP_jones[self.n_bands -(b+1)] is not None):
                           for chan in (a.BP_jones[self.n_bands -(b+1)]):
                              all_BP_xx.append(chan[pol_index]) 
                    else:
                        if(a.BP_jones[self.n_bands-(band_number-128)] is not None):
                            for chan in (a.BP_jones[self.n_bands-(band_number-128)]):
                                all_BP_xx.append(chan[pol_index])


        return all_BP_xx

    def BP_jones_amps(self,antenna=0):
        bp_jones_amps = []
        a = self.antennas[antenna]
        for B in a.BP_jones:
            bp_band_amps = []
            if(B!=None):
                for chan in B:
                    bp_band_amps.append(abs(chan))
            bp_jones_amps.append(bp_band_amps)
            
        return bp_jones_amps

    def all_BP_weights(self,cent_chan=12):
        all_weights = []
        for w in self.BP_weights:
            all_weights.append(w)

        return np.array(np.ravel(all_weights))    
    
    def single_fullband_jones(self,antenna=0,reverse_bands=False,correct_gain_jump=True,conjugate_jones=True,pol='xx',cent_chan=12,return_amps=False):

        all_single_jones = None
        a = self.antennas[antenna]
        band_0 = cent_chan - 12
        if(a!=None):
            for b,B in enumerate(a.BP_jones):
                band_number = b + cent_chan - 12
                if(band_number > 128):
                    if(band_0 > 128):
                        #entire band is reversed
                        if(a.Single_jones[self.n_bands -(b+1)] is not None):
                             if(all_single_jones is not None):
                                 all_single_jones = np.concatenate((all_single_jones,a.Single_jones[self.n_bands -(b+1)]))
                             else:
                                 all_single_jones = a.Single_jones[self.n_bands -(b+1)]
                        else:
                            if(a.Single_jones[self.n_bands-(band_number-128)] is not None):
                                if(all_single_jones is not None):
                                    all_single_jones = np.concatenate((all_single_jones,a.Single_jones[self.n_bands-(band_number-128)]))
                                else:
                                    all_single_jones = a.Single_jones[self.n_bands -(b+1)]
                else:
                    if(a.Single_jones[b] is not None):
                        if(all_single_jones is not None):
                            all_single_jones = np.concatenate((all_single_jones,a.Single_jones[b]))
                        else:
                            all_single_jones = a.Single_jones[b]

        return all_single_jones

    def all_fullband_jones(self,cent_chan=12):
        for i,a in enumerate(self.antennas):
            if(a.BP_jones[0]):
                a.full_band_jones = self.single_fullband_jones(antenna=i,cent_chan=cent_chan)

    def load_all_cals(self,path,metafile,n_bands=24,raw=True):

        # Check metafile exists
        if not os.path.exists(metafile):
            raise IOError(metafile + ' not found')

        self.load_metadata(metafile)

        self.load_all_BP_jones(path=path,raw=raw)
        self.load_all_DI_jones(path=path)

        self.form_single_jones()
        self.all_fullband_jones(cent_chan=self.cent_chan)
        
    def write_UCLA_cal(self,filename=None,reverse_bands=True):

        if(self.metafile is None):
            print("Error: Use load_metadata to define metafits file")
        else:    
            meta_file = fits.open(self.metafile)
        antenna_n = meta_file[1].data['Antenna']
        tile_n = meta_file[1].data['Tile']
        obsid = meta_file[0].header['GPSTIME']

        if(filename is None):
            out_file = open('RTS_%s_ucla.txt' % obsid,'w+')
        else:
            out_file = open(filename,'w+')

        # convert Cotter/FHD antenna number into rts antenna number
        ant2rts = dict(zip(antenna_n[::2],range(self.n_ants)))
        ant2tile = dict(zip(antenna_n[::2],tile_n[::2]))

        out_file.write("#Program of origin: RTS\n")
        out_file.write("# Convention: Divide uncalibrated data by these gains to obtain calibrated data.\n")

        pols = [['EE','EN'],['NE','NN']]

    #    for ai,a in enumerate(self.antennas):
        for p1 in [0,1]:
            for p2 in [0,1]:
                n_bands = 24
                n_chan = 16
                if(reverse_bands is True):
                    bands = reversed(range(n_bands))
                else:
                    bands = range(n_bands)
                for i_n,i in enumerate(bands):
                    for j in range(n_chan):
                        for k in range(self.n_ants):
                           if(ant2rts[k] not in self.flagged_antennas):
                               flag = 0
                               single_jones = self.antennas[ant2rts[k]].Single_jones[i][j]
                           else:
                               flag = 1
                               single_jones = np.zeros((2,2))
                               
                           out_file.write("%s %s %6.3f %s %.5f %f %f %d\n" % (ant2tile[k],k,self.freqs[i_n*n_chan+ j]/1.0e6,pols[p1][p2],self.mjd,np.real(single_jones[p1,p2]),np.imag(single_jones[p1,p2]),flag))

        out_file.close()

class rts_antenna():
    def __init__(self,ant_i,n_bands=24):
        self.ant_n = ant_i
        self.n_bands = n_bands
        self.BP_jones = [None] * n_bands
        self.DI_jones = [None] * n_bands
        self.Single_jones = [None] * n_bands
        self.full_band_jones = None
        self.DI_jones_ref = [None] * n_bands

    def average_BP_jones(self,fscrunch=2):
        for idx,B in enumerate(self.BP_jones):
            if(B!=None):
                avg_jones = [(B[2*i] + B[2*i+1])/2.0 for i in range(len(B)/2)]
                avg_jones[len(avg_jones)/2] *= 2.0
                self.BP_jones[idx] = avg_jones

    def form_single_jones(self,invert_order=False,flagged_channels=(0,1,16,30,31)):
        channels = np.arange(32)
        keep_channels = np.delete(channels,flagged_channels)
        
        for idx,[B,DI] in enumerate(zip(self.BP_jones,self.DI_jones)):
            if(B is not None):
                self.Single_jones[idx] = np.matmul(DI,B)
                if(invert_order):
                    self.Single_jones[idx] = [s * (np.matrix(self.DI_jones_ref[idx])).I for s in self.Single_jones[idx]]
                else:
                    invert_S = np.zeros_like(self.Single_jones[idx])
                    invert_S[keep_channels] = np.linalg.inv((self.Single_jones[idx])[keep_channels])
                    self.Single_jones[idx] = np.matmul(self.DI_jones_ref[idx],invert_S)                              

 
                    
def write_BP_files(raw_cal,fit_cal,filename='test',is_80khz=False,flagged_channels=[0,1,16,30,31]):
    from cmath import phase
    n_bands = 24
    n_tiles = 128
    if(is_80khz):
        n_chan = 16
        ch_width = 0.08
    else:
        n_chan = 32
        ch_width = 0.04
    bp_freqs = ""
    first_ch = True
    for ch in range(n_chan):
        if(ch not in flagged_channels):
            if(first_ch):
                bp_freqs = bp_freqs + "%7.6f" % (ch * ch_width)
                first_ch = False
            else:
                bp_freqs = bp_freqs + ", %7.6f" % (ch * ch_width)
                
    bp_freqs = bp_freqs + "\n"
        
    for n in range(raw_cal.n_bands):
        band_file = 'Bandpass' + filename + '%03d.dat' % (n+1)
        fp = open(band_file,'w+')
        fp.write(bp_freqs)
        for i in range(n_tiles):
            if(raw_cal.antennas[i] is not None and i not in raw_cal.flagged_antennas):
                if(raw_cal.antennas[i].BP_jones[n] is not None):
                    for pol1 in [0,1]:
                        for pol2 in [0,1]:
                            fp.write('%d' % (i+1))
                            for ch_n,ch in enumerate(raw_cal.antennas[i].BP_jones[n]):
                                if(ch_n not in flagged_channels):
                                    fp.write(', %f, %f' % (abs(ch[pol1,pol2]), phase(ch[pol1,pol2])))
                            fp.write('\n')
                            fp.write('%d' % (i+1))
                            for ch_n,ch in enumerate(fit_cal.antennas[i].BP_jones[n]):
                                if(ch_n not in flagged_channels):
                                    fp.write(', %f, %f' % (abs(ch[pol1,pol2]), phase(ch[pol1,pol2])))
                            fp.write('\n')
                    
        fp.close()

def write_DI_files(rts_cal,filename='test'):
    #n_bands = 24
    inv0_entry = "+1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0\n"
    for n in range(rts_cal.n_bands):
       band_file = 'DI_Jones' + filename + '%03d.dat' % (n+1)
       fp = open(band_file,'w+')
       fp.write("1.0\n")
       #       fp.write("%s" % inv0_entry)
       fp.write('%7.6f, %7.6f, %7.6f, %7.6f, %7.6f, %7.6f, %7.6f, %7.6f\n' % (np.real(rts_cal.antennas[0].DI_jones_ref[n][0,0]), np.imag(rts_cal.antennas[0].DI_jones_ref[n][0,0]),np.real(rts_cal.antennas[0].DI_jones_ref[n][0,1]), np.imag(rts_cal.antennas[0].DI_jones_ref[n][0,1]),np.real(rts_cal.antennas[0].DI_jones_ref[n][1,0]), np.imag(rts_cal.antennas[0].DI_jones_ref[n][1,0]),np.real(rts_cal.antennas[0].DI_jones_ref[n][1,1]), np.imag(rts_cal.antennas[0].DI_jones[n][1,1])))   
       
       for i,a in enumerate(rts_cal.antennas):
           if(a.DI_jones[n] is None or i in rts_cal.flagged_antennas):
             fp.write("+0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0\n") 
           else:
               fp.write('%7.6f, %7.6f, %7.6f, %7.6f, %7.6f, %7.6f, %7.6f, %7.6f\n' % (np.real(a.DI_jones[n][0,0]), np.imag(a.DI_jones[n][0,0]),np.real(a.DI_jones[n][0,1]), np.imag(a.DI_jones[n][0,1]),np.real(a.DI_jones[n][1,0]), np.imag(a.DI_jones[n][1,0]),np.real(a.DI_jones[n][1,1]), np.imag(a.DI_jones[n][1,1])))   
    

                
            
        
        
            
        
        


            
