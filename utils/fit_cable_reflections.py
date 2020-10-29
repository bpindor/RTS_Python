import numpy as np
from numpy import polyfit,poly1d, pi, arange, zeros,exp, pi, array, std, mean, random,real,imag, shape
import matplotlib.pyplot as plt
from scipy import optimize

# pads missing fine channels with zeros/gaussian noise for FFTs

def pad_edge_channels(bp_data,zero_padded=False,n_bands=24):

    sigma = np.std(bp_data)
    avg = np.mean(bp_data)

    flagged_channels = [0,1,16,30,31]

    chan_count = 0
    data_count = 0
    padded_bp = np.zeros(n_bands * 32)
    for i in range(n_bands):
        for j in range(32):
            if(j in flagged_channels):
                if(zero_padded):
                    padded_bp[data_count] = avg 
                else:
                    padded_bp[data_count] = avg + sigma * np.random.rand()
            else:
                padded_bp[data_count] = bp_data[chan_count]
                chan_count +=1
            data_count += 1

    return padded_bp

def clip_edge_channels(bp_data,clip_width=1,n_bands=24):

    chan_count = 0
    data_count = 0
    clipped_bp = np.zeros(n_bands * (27-2*clip_width),dtype=complex)
    clipped_x = np.zeros(n_bands * (27-2*clip_width))
    
    channels = np.arange(32)
    channels = np.delete(channels,16)
    channels = channels[2+clip_width:-(2+clip_width)]

    for i in range(n_bands):
        for j in range(32):
            if(j in channels):
                clipped_bp[data_count] = bp_data[chan_count]
                clipped_x[data_count] = chan_count
                data_count +=1
            chan_count += 1

    return clipped_bp, clipped_x

def bp_xvals(n_bands=24):
    # returns list of unflagged channel numbers
    flagged_channels = [0,1,16,30,31]
    chan_count = 0
    xvals = []
    for i in range(n_bands):
        for j in range(32):
            if(j not in flagged_channels):
                xvals.append(chan_count)
            chan_count +=1 
        
    return np.array(xvals)

def fit_reflection(data):

    sigma = np.std(data)
    avg = np.mean(data)
    phase = np.pi/2.0

    padded = pad_edge_channels(data)

    freq = 1.0 / (np.fft.fftfreq(len(padded))[np.argmax(abs(np.fft.fft(padded))[1:len(padded)/2])+1])

#    print freq
#    freq = 20.0
    #freq = 31.0
    freq = 24.0

    clipped, clipped_x = clip_edge_channels(padded,clip_width=4)

    sigma = np.std(clipped)
    avg = np.mean(clipped)

    #tX = data
    #Tx = arange(len(data))
    tX = clipped
    Tx = clipped_x

    p0 = [sigma,freq,phase,avg]
#    fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/p[1]*x+p[2]) + p[3] # Target function
    fitfunc = lambda p, x: (p[0]*sigma)*np.cos(2*np.pi/(p[1]*freq)*x+(p[2]*phase)) + (p[3]*avg)
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
    
#    print p0

    initial_vals = [1.0,1.0,0.0,1.0]

    p1, success = optimize.leastsq(errfunc, initial_vals, args=(Tx, tX))

#    print p1

    model_vals = [i * p for i,p in zip(p0,p1)]

#    print model_vals

    time = np.linspace(Tx.min(), Tx.max(), 1000)
#    plt.plot(Tx, tX, "ro", time, fitfunc(p1, time), "r-")
#    plt.plot(Tx, tX, "ro", time, fitfunc(p1, time), "r-")
#    plt.plot(time, fitfunc(p0, time), "b-")

    model_points = fitfunc(p1,bp_xvals())
    model_points = fitfunc(p0,bp_xvals())

    return model_points, model_vals

def fit_polynomial_bandpass(bp_data,n_bands=24,order=3):

    if(len(np.shape(bp_data)) == 1):

        xvals = np.arange(len(bp_data))
        
        clipped, clipped_x = clip_edge_channels(bp_data,clip_width=0,n_bands=n_bands)
    
        z_real = np.polyfit(clipped_x,np.real(clipped),order)
        z_imag = np.polyfit(clipped_x,np.imag(clipped),order)

        p_real = np.poly1d(z_real)
        p_imag = np.poly1d(z_imag)

        return p_real(xvals), p_imag(xvals)
    
    else:
        if(len(np.shape(bp_data)) == 3 and np.shape(bp_data)[0] == 2 and np.shape(bp_data)[1] == 2):
            xvals = np.arange(np.shape(bp_data)[2])
            models_out = [[[],[]],[[],[]]]
            for x in range(2):
                for y in range(2):
                    clipped, clipped_x = clip_edge_channels(bp_data[x][y],clip_width=0,n_bands=n_bands)  
                    z_real = np.polyfit(clipped_x,np.real(clipped),order)
                    z_imag = np.polyfit(clipped_x,np.imag(clipped),order)

                    p_real = np.poly1d(z_real)
                    p_imag = np.poly1d(z_imag)   
                    
                    models_out[x][y] = p_real(xvals) + 1.0j*p_imag(xvals) 

            return models_out

        else:
            print('Bandpass data of unsupported dimensions: %s' % np.shape(bp_data))
            exit(1)

            
def fit_auto_bandpass(data):

    sigma = np.std(data)
    avg = np.mean(data)
    phase = np.pi/2.0

#    data = pad_edge_channels(data)

    freq = 1.0 / (np.fft.fftfreq(len(data))[np.argmax(abs(np.fft.fft(data))[1:len(data)/2])+1])

    # Dont want to fit to channels with low pfb gain, so clip those
    
    clipped, clipped_x = clip_edge_channels(data,clip_width=4)

    sigma = np.std(clipped)
    avg = np.mean(clipped)

    #tX = data
    #Tx = arange(len(data))
    tX = clipped
    Tx = clipped_x

    p0 = [sigma,freq,phase,avg]
#    fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/p[1]*x+p[2]) + p[3] # Target function
    fitfunc = lambda p, x: (p[0]*sigma)*np.cos(2*np.pi/(p[1]*freq)*x+(p[2]*phase)) + (p[3]*avg)
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
    
#    print p0

    initial_vals = [1.0,1.0,0.0,1.0]

    p1, success = optimize.leastsq(errfunc, initial_vals, args=(Tx, tX))

#    print p1

    model_vals = [i * p for i,p in zip(p0,p1)]

#    print model_vals

    time = np.linspace(Tx.min(), Tx.max(), 1000)
#    plt.plot(Tx, tX, "ro", time, fitfunc(p1, time), "r-")
#    plt.plot(Tx, tX, "ro", time, fitfunc(p1, time), "r-")
#    plt.plot(time, fitfunc(p0, time), "b-")

    model_points = fitfunc(p1,bp_xvals())
    model_points = fitfunc(p0,bp_xvals())

    return model_points, model_vals
    
def fit_polynomial_and_reflection_bandpass(bp_data,n_bands=24,order=3,return_seed=False):

    # Dont want to fit to channels with low pfb gain, so clip those
    
    clipped, clipped_x = clip_edge_channels(bp_data,clip_width=4)

    tX = np.real(clipped)
    Tx = clipped_x

    # How to seed wavelength, phase??

    freq = 300.0
    phase = np.pi/2.0

    sigma = np.std(tX)
    avg = np.mean(tX)

    # Scaling factors so fit parameters ~ 1
    
    p0 = [sigma,freq,phase,avg]
    initial_vals = [1.0,1.0,0.0,1.0]

    print(p0)

    fitfunc = lambda p, x: (p[0]*sigma)*np.cos(2*np.pi/(p[1]*freq)*x+(p[2]*phase)) + (p[3]*avg)
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
    p1, success = optimize.leastsq(errfunc, initial_vals, args=(Tx, tX))

#    print p1

    model_vals = [i * p for i,p in zip(p0,p1)]

#    print model_vals

    time = np.linspace(Tx.min(), Tx.max(), 1000)

    model_points = fitfunc(p1,bp_xvals())
    if(return_seed):
        model_points = fitfunc(initial_vals,bp_xvals())
    else:
        model_points = fitfunc(p1,bp_xvals())

    return model_points, model_vals


    
