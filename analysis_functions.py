# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 15:35:08 2017
@author: dplatzer

Some of the treatment functions by Margherita, almost unchanged. Contains:
find_local_maxima
tof2Eev
jacobian_transform
smooth
FFT
wrap2pmpi
butter_bandpass
find_2w
"""

from scipy.signal import gaussian, hamming, butter, freqz
from glob import glob
from numpy import loadtxt, arange, array, where, zeros, convolve
from numpy import savetxt, int32, append, linspace, argmin
import numpy as np
from scipy import optimize as opt
from scipy import interpolate
import glob_var as cts


def find_local_maxima(d, thry, thrxmin, thrxmax, smooth1=5, smooth2=20, mindt=1):
	"""Find local maxima of a series. Each peak is found by adaptive
    thresholding.
    Return the indices of the local maxima and the value at that index.
    The threshold is adaptive.
    Parameters:
     * smooth1 serves to smooth the series before finding the peaks. A value
    around 2-5 should be ok.
     * smooth2 serves for the adaptive thresholding: we consider a peak all
    what is above the data smoothed by smooth2. Choose something a bit smaller
    than the duration of a pulse
     * mindt serves to count as a single peak the ones which are
    separated by less than mindt. For a smooth series where the thresholding
    works well, the default value is ok. If the peaks are not well behaved and
    have a bit of digitalization noise, set mindt to something close to the
    width of a peak and it should solve the problem.


       Example:
           a,b=ts.find_local_maxima( dd, 0)
           plot(dd,gpdata(a,b))

    NOTE: this could be enhanced by multiplying the adaptive threshold by a
    value bigger than 1. It would remove false positives. Think about it"""
# Threshold the series: we threshold a smoothed series to find the peaks
# approximately by adaptive thresholding
	ds=convolve(d,gaussian(2**12,smooth1)/gaussian(2**12,smooth1).sum(),mode='same')
#	ds=d
	thresh=convolve(d,gaussian(2**12,smooth2)/gaussian(2**12,smooth2).sum(),mode='same')

	dt = array(where(ds[0:len(thresh)]> thresh))

	dt = dt[0,:]

    # Compute the derivative of the thresholded series: where the
    # derivative  is big means that we jumped to a new peak
	w = dt[1:]-dt[0:-1]

    # Label the peaks
	peaknumber = 0
	wp = zeros(w.shape, int32)
	for i in range(0, len(w)):
		if w[i] <= mindt:
			wp[i] = peaknumber
		else:
			wp[i] = peaknumber
			peaknumber = peaknumber + 1

    # Find the index of the maximum inside each peak
	peaklocations = zeros(wp.max(), int32)
 #   peaklocations2 = zeros(wp.max(),int32)
	for i in range(0,wp.max()):
		aux = dt[where(wp == i)]
#	        peaklocations[i] =  d[list(aux)].argmax()+aux[0]
		peaklocations[i] =  ds[list(aux)].argmax()+aux[0]
#
	dpeak=[]
	ipeak=[]

#

	for i in range(0,len(peaklocations)):
		if ds[peaklocations[i]] > thry and peaklocations[i]*1e-9 >= thrxmin and peaklocations[i]*1e-9 <= thrxmax: #modified by dom
			dpeak.append(ds[peaklocations[i]])
			ipeak.append(peaklocations[i])

	dpeak=array(dpeak)
	ipeak=array(ipeak)

	#plot(ipeak,dpeak,'ro')
	return ipeak,dpeak,ds

def tof2EeV(tof,a,t0,c):
    """ Function to convert time of flight in Energy [eV]"""
    return a/(tof-t0)**2 + c

def jacobian_transform(tof, counts):

    a_fit = cts.afit
    t0_fit = cts.t0fit
    c_fit = cts.cfit
    nu = cts.cur_nu
    dE = cts.dE
    elow = cts.elow
    ehigh = cts.ehigh

    qq = a_fit/(tof-t0_fit)**2 + c_fit
    EeV = cts.HEV*nu*qq
    qqderiv = 2*a_fit/(tof-t0_fit)**3
    Ederiv = abs(cts.HEV*nu*qqderiv)
    #print(Ederiv)
#    a_fit,t0_fit,c_fit,nu = np.loadtxt(fdir+'calib.txt')
#    EeV = a_fit/(tof-t0_fit)**2 + c_fit
#    Ederiv = 2*a_fit/(tof-t0_fit)**3


    eevlin = np.arange(elow,ehigh,step=dE)


    dshape = counts.shape
    #print(dshape)

    if len(dshape) > 1:
        counts = counts[:,::-1]
        dshape = counts.shape
        countsnew = np.zeros((dshape[0],dshape[1]))

        for i in range(dshape[0]):
            countsnew[i,:]=counts[i,:]/Ederiv

        f = interpolate.interp1d(EeV,countsnew)
        countsnew = f(eevlin)
        print("blop")
    else:
        counts = counts[::-1]
        countsj = counts/Ederiv
        countsnew = np.interp(eevlin,EeV,countsj)
    #print(countsnew)

    return eevlin, countsnew

def smooth(data,sm=2):
    """It smooths data by convolvolving it with a gaussian function of width sm"""
    nn = len(data)
    nexp = [n for n in range(5,15) if 2**n < nn][-1]

    ds=np.convolve(data,gaussian(2**nexp,sm)/gaussian(2**nexp,sm).sum(),mode='same')
    return ds
#def tof_to_energy(fitfn,):
def FFT(x_data, dt):

     nu = cts.cur_nu
     zero_order = cts.FT_zero_order
     window = cts.FT_window
     padding = cts.FT_padding
     npad = cts.FT_npad

     nu = nu*1e-15

     ##############
     N = len(x_data) # number of vector points
#
     n = npad  # number of points for padding
     pad = np.zeros(n)
     # use only even number of points for the signal
     if N % 2 == 0:
         N = N
     else:
         N = N-1
         x_data =  x_data[1:]


     # window for limiting spectral leakage and avoid sidelobes in fft amplitude
     Hm = hamming(N)


     signal =  x_data
     if zero_order == False:
         signal = signal - signal.mean()

     if window == True:
         #print('Using Hamming window')
         signal = signal*Hm
     else:
         signal = signal

     if padding == True:
         signal = signal
         ###### NAL 000 SIG ###
         pad[:N//2] = signal[N//2:]
         pad[-N//2:] = signal[:N//2]
         ###### 000 SIGNAL 000 ###
         	# a = n - N
         	# pad[a//2:-a//2] = F1
         signal=pad
         fourier = np.fft.fft(signal)
         fourier = fourier[:n//2]

         freq = np.fft.fftfreq(n, d=dt)
         freqpos = freq[:n//2]
         freqnorm = freqpos/nu
     else:
         signal = signal

         fourier = np.fft.fft(signal)
         fourier = fourier[:N//2]

         freq = np.fft.fftfreq(N, d=dt)
         freqpos = freq[:N//2]
         freqnorm = freqpos/nu

     ampl= abs(fourier)


     phicorr = -np.pi*freqpos*dt

     ang = np.angle(fourier) + phicorr
     ang=wrap2pmpi(ang)


     return freqnorm,ampl,ang

def wrap2pmpi(phasedata):
    """It wraps phase from -pi, pi"""
    return (phasedata + np.pi) % (2*np.pi) -np.pi

def butter_bandpass(lowcut, highcut, fs, n, order=5):
    """It builds a butter bandpass filter.
            Parameters: lowcut = low frequency
                        highcut = high frequency
                        fs = sampling frequency
                        n = number of points
                        order = order of butter filter

            It returns: w = the normalized frequencies at which h is computed
                        h = the frequency response """
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='bandpass')
    w, h = freqz(b,a,worN = n)
    return w, h



def find_2w(x_data, dt):

    nu = cts.cur_nu
    average = cts.two_w_average
    bfilter = cts.two_w_bfilter
    integral = cts.two_w_integral
    phi_offset = cts.two_w_phioffset
   ###### Sampling frequency and bandwidth for the butter bandpass filter ########
    fs = 1/dt/nu*1e15
    BW = 1.2
     ##############
#
    N = len(x_data)

#### real frequency step of fft(signal)######
    dfreal = 1/(dt*N)/nu*1e15

    freqnorm,ampl,ang=FFT(x_data,dt)
       ####### find index corresponding to the 2-omega peak #########
    pk_indeces = np.where(abs(freqnorm - 2) < 0.5)[0]

    peak_index = np.argmax(ampl[pk_indeces])
    peak_index_good = pk_indeces[peak_index]
    fpeak = freqnorm[peak_index_good]
    #ifmin=ifmax=0
     ########### filtered fft ############
    if bfilter == True and average == False:
        #print('Using bandpass filter around 2omega')
        lowcut = fpeak - BW/2.
        highcut = fpeak + BW/2.
         #i_lc = np.argmin(abs(freqnorm - lowcut))
         #i_hc = np.argmin(abs(freqnorm - highcut))
        n=len(freqnorm)
        w, h = butter_bandpass(lowcut,highcut, fs, n, 5)
        ampl_flt = ampl*h
#             power_spectrum = abs(fourier_flt)**2

        ang_flt = ang + phi_offset
#             ang_flt = np.unwrap(ang_flt)
        peak_phase = ang_flt[peak_index_good]
        peak = ampl_flt[peak_index_good]

    elif bfilter == False and average == True:
        #print('Averaging phase over the real frequency step')
#     ########### average phase over dfreal #########
        fmin = fpeak - dfreal/2
        fmax = fpeak + dfreal/2
        ifmin = np.argmin(abs(freqnorm-fmin))
        ifmax = np.argmin(abs(freqnorm-fmax))

        if(ifmin<ifmax):
           peak = ampl[ifmin:ifmax].mean()
           peak_phase = np.mean(np.unwrap(ang[ifmin:ifmax])) + phi_offset
        else:
           peak = ampl[ifmin]
           peak_phase = 	ang[ifmin] + phi_offset
    elif bfilter == True and average == True:
        #print('Using bandpass filter and averaging over the real frequency step')
        lowcut = fpeak - BW/2.
        highcut = fpeak + BW/2.
        n=len(freqnorm)
        w, h = butter_bandpass(lowcut,highcut, fs, n, 5)
        ampl_flt = ampl*h
        ang_flt = ang + phi_offset
        fmin = fpeak - dfreal/2
        fmax = fpeak + dfreal/2
        ifmin = np.argmin(abs(freqnorm-fmin))
        ifmax = np.argmin(abs(freqnorm-fmax))
        if(ifmin<ifmax):
           peak = ampl[ifmin:ifmax].mean()
           peak_phase = np.mean(np.unwrap(ang[ifmin:ifmax])) + phi_offset
        else:
           peak = ampl[ifmin]
           peak_phase = 	ang[ifmin] + phi_offset
    else:
        #print('Taking the peak phase')
        peak = ampl[peak_index_good]
        peak_phase = ang[peak_index_good] + phi_offset

    if integral == True and average == True:
        print('Error! either average or integral, not both!')
    elif integral == True and bfilter == True:
        fourier_flt = ampl_flt*np.exp(1j*ang_flt)
        fourier_flt_int = np.trapz(fourier_flt,freqnorm)
        peak_phase = np.angle(fourier_flt_int) + phi_offset
        peak = abs(fourier_flt_int)

    #print(ifmin,ifmax)
    #print(fpeak,peak,peak_phase)

    return fpeak,peak,peak_phase





