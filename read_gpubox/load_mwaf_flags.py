from astropy.io import fits 
import sys,os
from numpy import shape, sum, zeros

def load_mwaf_flags(mwaffile):
    
    hdulist = fits.open(mwaffile)
    nchan=hdulist[0].header['NCHANS']
    nant =hdulist[0].header['NANTENNA']
    nbl  = nant*(nant+1)/2
    
    flags = hdulist[1].data['FLAGS']

#There can be missing time samples so
#ntime=hdulist[0].header['NSCANS'], instead
    ntime = shape(flags)[0] / nbl

    all_flags = flags.reshape(ntime,nbl,nchan)

    return all_flags


