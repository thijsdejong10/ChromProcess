import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
import matplotlib.pyplot as plt
def _1gaussian(x, amplitude, centre, sigma):
    '''
    A single gaussian function

    Parameters
    ----------
    x: array
        x axis data
    amplitude1: float
        amplitudelitude of the function
    centre1: float
        centretre of the function (mean)
    sigma1: float
        width of the function (standard deviation)
    Returns
    -------
    function: numpy array
        y values for the function
    '''
    return amplitude*(1/(sigma*(np.sqrt(2*np.pi))))*(np.exp(-((x-centre)**2)/((2*sigma)**2)))

def _2gaussian(x, amplitude1, centre1, sigma1, amplitude2, centre2, sigma2):
    '''
    A double gaussian function
    Parameters
    ----------
    x: array
        x axis data
    amplituden: float
        amplitudelitude of a component gaussian function
    centren: float
        centretre of a component gaussian function (mean)
    sigman: float
        width of a component gaussian function (standard deviation)
    Returns
    -------
    function: numpy array
        y values for the function
    '''
    return _1gaussian(x, amplitude1, centre1, sigma1) + _1gaussian(x, amplitude2, centre2, sigma2)

def _3gaussian(x, amplitude1, centre1, sigma1,amplitude2, centre2, sigma2, amplitude3, centre3, sigma3):
    return _1gaussian(x, amplitude1, centre1, sigma1) + _2gaussian(x, amplitude2, centre2, sigma2, amplitude3, centre3, sigma3)

def fit_gaussian_peaks(
                    time, 
                    signal,
                    initial_guess,
                    boundaries,
                    num_peaks,
                    ):
    '''
    TODO: This kind of function could be useful, but a better adapted function
        for peak deconvolution should be written. The function could take
        similar concepts to this one, but with a different concept for its
        implementation.

    Fitting sums of gaussian peaks to data using supplied peak indices.

    Parameters
    ----------
    time: array
        time values
    sig: array
        signal to be fit to
    peaks: list of peak positions
        list peak positions in the time

    initial_guess: list
        Initial guess for the peak amplitudelitude, position and width
        e.g. see _1gaussian() function arguments.

    lowerbounds: list
        Lower bounds for the peak amplitudelitude, position and width
        e.g. see _1gaussian() function arguments.
    upperbounds: list
        Upper bounds for the peak amplitudelitude, position and width
        e.g. see _1gaussian() function arguments.

    Returns
    -------
    popt: ndarray
        list of fitted values [[amplitudelitude, centretre, width],]
    pcov: array
        correlation matrix
    '''

    if num_peaks == 1:
        popt, pcov = curve_fit(_1gaussian, time, signal, p0=initial_guess, bounds = boundaries)
    elif num_peaks == 2:
        popt, pcov = curve_fit(_2gaussian, time, signal, p0=initial_guess, bounds = boundaries)
    elif num_peaks == 3:
        popt, pcov = curve_fit(_3gaussian, time, signal, p0=initial_guess, bounds = boundaries)#,method='trf')
    else:
        print('Error: number of peaks to large')
        popt, pcov = [0,0,0],0

    return popt, pcov

def deconvolute_region(chromatogram, region, peak_diff, num_peaks = 1):

    '''
    TODO: Combine the ideas in this function with fit_gaussian_peaks()

    chromatogram: ChromProcess Chromatogram object
        Chromatogram
    region: list
        region of chromatogram under operation [lower bound, upper bound]
    '''

    upper = region[1]
    lower = region[0]

    inds = np.where((chromatogram.time > lower)&(chromatogram.time < upper))[0]

    time = chromatogram.time[inds]
    signal = chromatogram.signal[inds]

    signal = signal-np.average(signal[-5:-1])

    for p in range(0,len(peaks)): # extend the boundary
        initial_guess[0] = np.amax(signal)
        initial_guess[1] = peaks[p]
        lowerbounds[1] = time[0]
        upperbounds[1] = time[-1]
        guess.extend(initial_guess)
        lbds.extend(lowerbounds)
        ubds.extend(upperbounds)

    #peak_list = np.array([*chromatogram.peaks])
    #peak_inds = np.where((peak_list > lower)&(peak_list < upper))[0]
#
    #peaks = peak_list[peak_inds]
#
    #while len(peaks) < num_peaks:
    #    peaks = np.append(peaks, np.average(peaks))
#
    #if len(peaks) > num_peaks:
    #    peaks = peaks[:num_peaks]


    return fit_gaussian_peaks(time, signal, peaks)

def deconvolute_peak(
                    chromatogram, 
                    initial_guess,
                    experiment_folder,
                    indices,
                    lower_bounds = [0, 0,  0.0], 
                    upper_bounds = [1e10, 1,  0.035], 
                    num_peaks = 2,
                    plotting = True,
                    ):

    '''
    TODO: this function is quite similar in scope to deconvolute_region().
    Refactor with the other two deconvolution macros.

    chromatogram: ChromProcess Chromatogram object
        Chromatogram
    region: list
        region of chromatogram under operation [lower bound, upper bound]
    '''
    from ChromProcess.Utils.signal_processing import signal_processing as sig
    time = chromatogram.time[indices[0]:indices[-1]]
    signal = chromatogram.signal[indices[0]:indices[-1]]
    smoothed_signal = sig.savitzky_golay(signal, 5, 3, deriv=0, rate=1)
    #lower_bounds = []
    #upper_bounds = []
    #for p in range(0,num_peaks): # extend the boundary
    #    standard_lower_bounds[0] = min(signal)*0.1
    #    standard_lower_bounds[1] = time[0]
    #    standard_upper_bounds[0] = max(signal)*3
    #    standard_upper_bounds[1] = time[-1]
    #    lower_bounds.extend(standard_lower_bounds)
    #    upper_bounds.extend(standard_upper_bounds)
    #
    boundaries = [lower_bounds,upper_bounds]

    popt, pcov  =  fit_gaussian_peaks(time, signal, initial_guess, boundaries, num_peaks)
    
    if plotting==True:
        fig, ax = plt.subplots(1,1)
        ax.plot(time,signal)
        final_curve = np.zeros(len(signal))
        for p in range(0,num_peaks):
            gauss = _1gaussian(time,popt[0+p*3],popt[1+p*3],popt[2+p*3])
            final_curve += gauss
            ax.plot(time,gauss)
        ax.plot(time,final_curve)
        fig.set_size_inches(18.5, 10.5)
        plt.tight_layout()

        plt.savefig(f"{experiment_folder}\\deconvolve_peaks\\11.70\\{chromatogram.filename}.png")
        plt.close()
    
    return popt, pcov 
