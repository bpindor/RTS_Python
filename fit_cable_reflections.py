import numpy as np
#from numpy import polyfit,poly1d, pi, arange, np.zeros,exp, pi, array, std, mean, random,real,imag, shape
import matplotlib.pyplot as plt
from scipy import optimize

# pads missing fine channels with np.zeros/gaussian noise for FFTs

def pad_edge_channels(bp_data,zero_padded=False):

    sigma = std(bp_data)
    avg = mean(bp_data)

    flagged_channels = [0,1,16,30,31]

    chan_count = 0
    data_count = 0
    padded_bp = np.zeros(24 * 32)
    for i in range(24):
        for j in range(32):
            if(j in flagged_channels):
                if(zero_padded):
                    padded_bp[data_count] = avg 
                else:
                    padded_bp[data_count] = avg + sigma * random.rand()
            else:
                padded_bp[data_count] = bp_data[chan_count]
                chan_count +=1
            data_count += 1

    return padded_bp

def clip_edge_channels(bp_data,clip_width=1,n_bands=24):

    chan_count = 0
    data_count = 0
    #clipped_bp = np.zeros(n_bands * (27-2*clip_width),dtype=complex)
    clipped_bp = np.zeros(n_bands * (27-2*clip_width))
    clipped_x = np.zeros(n_bands * (27-2*clip_width))
    
    channels = range(32)
    channels.pop(16)
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

#    p0 = [sigma,freq,phase,avg]
#    fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/p[1]*x+p[2]) + p[3] # Target function
#    fitfunc = lambda p, x: (p[0]*sigma)*np.cos(2*np.pi/(p[1]*freq)*x+(p[2]*phase)) + (p[3]*avg)
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

def fit_bandpass_and_reflection(bp_data,n_bands=24,order=7):

    if(len(np.shape(bp_data)) == 1):

        xvals = np.arange(len(bp_data))
        
        clipped, clipped_x = clip_edge_channels(bp_data,clip_width=4,n_bands=n_bands)
    
        z_real = np.polyfit(clipped_x,np.real(clipped),3)
        z_imag = np.polyfit(clipped_x,np.imag(clipped),3)

        p_real = np.poly1d(z_real)
        p_imag = np.poly1d(z_imag)

        return p_real(xvals), p_imag(xvals)
    
    else:
        if(len(np.shape(bp_data)) == 3 and np.shape(bp_data)[0] == 2 and np.shape(bp_data)[1] == 2):
            xvals = np.arange(np.shape(bp_data)[2])
            models_out = [[[],[]],[[],[]]]
            for x in range(2):
                for y in range(2):
                    clipped, clipped_x = clip_edge_channels(bp_data[x][y],clip_width=0)  
                    z_real = np.polyfit(clipped_x,np.real(clipped),3)
                    z_imag = np.polyfit(clipped_x,np.imag(clipped),3)

                    p_real = np.poly1d(z_real)
                    p_imag = np.poly1d(z_imag)   
                    
                    models_out[x][y] = p_real(xvals) + 1.0j*p_imag(xvals) 

            return models_out

        else:
            print 'Bandpass data of unsupported dimensions: %s' % np.shape(bp_data)
            exit(1)

def residuals(params,*args):


    data = args[0]
    x_vals = args[1]
    scale = args[2]
    freq = args[3]
    phase = args[4]
    y_int = args[5]
    slope = args[6]
    quad = args[7]
    cubic = args[8]
    
    fitfunc = lambda p, x: (p[0]*scale)*np.cos(2*np.pi/(p[1]*freq)*x+(p[2]*phase)) + (p[3]*y_int) + (p[4]*slope)*x + (p[5]*quad)*x*x + (p[6]*cubic)*x*x*x + (p[7]-1.0)*x*x*x*x + (p[8]-1.0)*x*x*x*x*x
    data = args[0]
    x_vals = args[1]
    
    model = fitfunc(params,x_vals)
    residuals = sum(pow((model-data),2.0))
    return residuals
            
def fit_auto_bandpass3(data,n_bands=24):

    
    clipped, clipped_x = clip_edge_channels(data,clip_width=4,n_bands=n_bands)

    # First fit low order polynomial
    
    z_real = np.polyfit(clipped_x,np.real(clipped),3)

    coarse_model = np.poly1d(z_real)

    coarse_clipped = coarse_model(clipped_x)

    # Divide out to estimate sin wave parameters

    coarse_ratio = clipped / coarse_clipped

    scale = max(coarse_ratio) - min(coarse_ratio)
    print scale

    scale *= np.mean(coarse_model)
    print scale

    #Still need to tune these
    
    phase = - np.pi / 2.0
    freq = 250
    
    ##########

    y_int = z_real[3]
    slope = z_real[2]
    quad  = z_real[1]
    cubic = z_real[0]

    min_info = (clipped,clipped_x,scale,freq,phase,y_int,slope,quad,cubic)

    initial_vals = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]

    res = optimize.minimize(residuals, initial_vals, args=min_info,method='Nelder-Mead',options={'disp': True, 'maxfev': 30000, 'maxiter': 10000})
#    res = optimize.minimize(residuals, initial_vals, args=(data),method='Nelder-Mead',options={'disp': True})

    print res.x
    print res.nit

    fitfunc = lambda p, x: (p[0]*scale)*np.cos(2*np.pi/(p[1]*freq)*x+(p[2]*phase)) + (p[3]*y_int) + (p[4]*slope)*x + (p[5]*quad)*x*x + (p[6]*cubic)*x*x*x + (p[7]-1.0)*x*x*x*x + (p[8]-1.0)*x*x*x*x*x

    model_points = fitfunc(res.x,np.arange(len(data)))
    initial_points = fitfunc(initial_vals,np.arange(len(data)))

    return model_points,res.x,initial_points

    
            
def fit_auto_bandpass2(data,n_bands=24):

    
    clipped, clipped_x = clip_edge_channels(data,clip_width=4,n_bands=n_bands)

    # First fit low order polynomial
    
    z_real = np.polyfit(clipped_x,np.real(clipped),3)

    coarse_model = np.poly1d(z_real)

    coarse_clipped = coarse_model(clipped_x)

    # Divide out to estimate sin wave parameters

    coarse_ratio = clipped / coarse_clipped

    scale = max(coarse_ratio) - min(coarse_ratio)
    print scale

    scale *= np.mean(coarse_model)
    print scale

    #Still need to tune these
    
    phase = - np.pi / 2.0
    freq = 250
    
    ##########

    y_int = z_real[3]
    slope = z_real[2]
    quad  = z_real[1]
    cubic = z_real[0]

    #initial_vals = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
    #p0 = [scale,freq,phase,y_int,slope,quad,cubic,1.0,1.0,1.0,1.0]
    initial_vals = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
    p0 = [scale,freq,phase,y_int,slope,quad,cubic,1.0,1.0]
    
    
    tX = clipped
    Tx = clipped_x

    print 'slope',slope

    #scale=0
    
    #fitfunc = lambda p, x: (p[0]*scale)*np.cos(2*np.pi/(p[1]*freq)*x+(p[2]*phase)) + (p[3]*y_int) + (p[4]*slope)*x + (p[5]*quad)*x*x + (p[6]*cubic)*x*x*x + (p[7]-1.0)*x*x*x*x + (p[8]-1.0)*x*x*x*x*x + (p[9]-1.0)*x*x*x*x*x*x + + (p[10]-1.0)*x*x*x*x*x*x*x
    fitfunc = lambda p, x: (p[0]*scale)*np.cos(2*np.pi/(p[1]*freq)*x+(p[2]*phase)) + (p[3]*y_int) + (p[4]*slope)*x + (p[5]*quad)*x*x + (p[6]*cubic)*x*x*x + (p[7]-1.0)*x*x*x*x + (p[8]-1.0)*x*x*x*x*x
    
#    fitfunc = lambda p, x: (p[0]*scale)*np.cos(2*np.pi/(p[1]*freq)*x+(p[2]*phase)) + (p[3]*y_int) + (p[4]*slope)*x + (p[5]-1.0)*x*x + (p[6]-1.0)*x*x*x + (p[7]-1.0)*x*x*x*x + (p[8]-1.0)*x*x*x*x*x + (p[9]-1.0)*x*x*x*x*x*x + + (p[10]-1.0)*x*x*x*x*x*x*x
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
    p1, success = optimize.leastsq(errfunc, initial_vals, args=(Tx, tX))
    full_output = optimize.leastsq(errfunc, initial_vals, args=(Tx, tX),full_output=True)
    

    print p0
    print p1
    print 'success', success
    
    model_vals = [i * p for i,p in zip(p0,p1)]

    print model_vals
    
    model_points = fitfunc(p1,np.arange(len(data)))
    initial_points = fitfunc(initial_vals,np.arange(len(data)))

    return model_points, model_vals, initial_points
    #return full_output, coarse_clipped

            
def fit_auto_bandpass(data,n_bands=24):

#    sigma = np.std(data)
#    avg = np.mean(data)
#    phase = np.pi/2.0
#    phase = 0.0

    sigma = 5000.0
    phase = - np.pi / 2.0
    freq = 25
    avg = 50000.0
    slope = (data[-1] - data[0]) / float(len(data))

#    print slope
    
#    data = pad_edge_channels(data)

#    freq = 1.0 / (np.fft.fftfreq(len(data))[np.argmax(abs(np.fft.fft(data))[1:len(data)/2])+1])

    freq = 250

    # Dont want to fit to channels with low pfb gain, so clip those
    
    clipped, clipped_x = clip_edge_channels(data,clip_width=4,n_bands=n_bands)

    sigma = np.std(clipped)
    avg = np.mean(clipped)
    slope = (clipped[-1] - clipped[0]) / float(len(data))
    y_int = clipped[0]

    #tX = data
    #Tx = arange(len(data))
    tX = clipped
    Tx = clipped_x

#    tX = np.array([d for d in data])
#    Tx = np.arange(len(data))

#    tX = sigma * np.cos(2*np.pi/(freq)*Tx+(phase)) + (avg)
    
#    p0 = [sigma,freq,phase,avg,slope,1.0,1.0,1.0,1.0]
#    p0 = [sigma,freq,phase,avg,slope,1.0,1.0,1.0,1.0]
    p0 = [sigma,freq,phase,y_int,slope,1.0,1.0,1.0,1.0,1.0,1.0]

#    fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/p[1]*x+p[2]) + p[3] # Target function
#    fitfunc = lambda p, x: (p[0]*sigma)*np.cos(2*np.pi/(p[1]*freq)*x+(p[2]*phase)) + (p[3]*avg) + (p[4]*slope)*x + (p[5]-1.0)*x*x + (p[6]-1.0)*x*x*x + (p[7]-1.0)*x*x*x*x + (p[8]-1.0)*x*x*x*x*x
    fitfunc = lambda p, x: (p[0]*sigma)*np.cos(2*np.pi/(p[1]*freq)*x+(p[2]*phase)) + (p[3]*y_int) + (p[4]*slope)*x + (p[5]-1.0)*x*x + (p[6]-1.0)*x*x*x + (p[7]-1.0)*x*x*x*x + (p[8]-1.0)*x*x*x*x*x + (p[9]-1.0)*x*x*x*x*x*x + + (p[10]-1.0)*x*x*x*x*x*x*x
    

#    p0 = [
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
#    print p0

#    initial_vals = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
    initial_vals = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]

#    initial_vals = (1.0,1.0,1.0,1.0)

#    p1, success = optimize.leastsq(errfunc, initial_vals, args=(Tx, tX))
    p1, success = optimize.leastsq(errfunc, initial_vals, args=(Tx, tX))

#    p1 = initial_vals
    
#    print p1

    model_vals = [i * p for i,p in zip(p0,p1)]

#    print model_vals

    time = np.linspace(Tx.min(), Tx.max(), 1000)
#    plt.plot(Tx, tX, "ro", time, fitfunc(p1, time), "r-")
#    plt.plot(Tx, tX, "ro", time, fitfunc(p1, time), "r-")
#    plt.plot(time, fitfunc(p0, time), "b-")

    #model_points = fitfunc(p1,bp_xvals(n_bands=n_bands))
    #model_points = fitfunc(p0,bp_xvals(n_bands=n_bands))
    model_points = fitfunc(p1,np.arange(len(data)))
    initial_points = fitfunc(initial_vals,np.arange(len(data)))
    
    return model_points, model_vals, initial_points
    
def fit_auto_bandpass_log_poly(data,n_bands=24):

    clipped, clipped_x = clip_edge_channels(data,clip_width=4,n_bands=n_bands)
    slope = (clipped[-1] - clipped[0]) / float(len(data))
    y_int = clipped[0]
    
    print slope
    
    tX = clipped
    Tx = clipped_x

    ltX = np.log10(clipped)
    lTx = np.log10(clipped_x)
    
#    p0 = [sigma,freq,phase,avg,slope,1.0,1.0,1.0,1.0]

    fitfunc = lambda p, x: p[0] + p[1] * x
    powerfunc = lambda p, x: p[0] * pow(x,p[1])
    
#    fitfunc = lambda p, x: (p[0]*sigma)*np.cos(2*np.pi/(p[1]*freq)*x+(p[2]*phase)) + (p[3]*avg) + (p[4]*slope)*x + (p[5]-1.0)*x*x + (p[6]-1.0)*x*x*x + (p[7]-1.0)*x*x*x*x + (p[8]-1.0)*x*x*x*x*x

    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
    

    initial_vals = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]

    initial_vals = [1.0, -1.0]
    
    #    p1, success = optimize.leastsq(errfunc, initial_vals, args=(Tx, tX))

    p1, success = optimize.leastsq(errfunc, initial_vals, args=(lTx, ltX))  

    model_vals = [i * p for i,p in zip(p0,p1)]


    time = np.linspace(Tx.min(), Tx.max(), 1000)

    model_points
    
    model_points = powerfunc(p1,np.arange(len(data)))
    initial_points = powerfun(initial_vals,np.arange(len(data)))
    
    return model_points, model_vals, initial_points

def get_pfb_from_bp(bp_data,n_bands=24):

    n_fine = 32
    pfb_shape = np.zeros(n_fine)

    for n in range(n_bands):
        pfb_shape += bp_data[n*n_fine:(n+1)*n_fine]

    pfb_shape /= (float) (n_bands)

    return pfb_shape

def apply_pfb_to_autos(auto_data,flags,pfb_model,n_bands=24):

    n_times = (auto_data.shape)[0]
    n_ants = (auto_data.shape)[2]
    n_fine = 32
    data_out = np.zeros(auto_data.shape)
    
    for t in range(n_times):
        for n in range(n_ants):
            if(flags[n] == 0):
                for b in range(n_bands):
                    data_out[t,b*n_fine:(b+1)*n_fine,n] = auto_data[t,b*n_fine:(b+1)*n_fine,n] / pfb_model

    return data_out


    
