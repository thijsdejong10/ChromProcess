import os
from statistics import mean
import matplotlib.pyplot as plt
from ChromProcess.Loading import chrom_from_csv
from ChromProcess.Loading import analysis_from_csv
from ChromProcess.Loading import conditions_from_csv
from ChromProcess.Loading.peak.peak_from_csv import peak_rt_from_file, peak_boundaries_from_file
experiment_number = 'FRN140'
experiment_folder = r"C:\Users\thijs\Documents\PhD\Data\FRN140"
from ChromProcess.Utils.peak_finding import find_peaks_scipy
from ChromProcess.Utils import indices_from_boundary, peak_indices_to_times
from ChromProcess.Utils.signal_processing.deconvolution import _1gaussian, _2gaussian, _3gaussian, deconvolute_peak, fit_gaussian_peaks
from ChromProcess.Processing import add_peaks_to_chromatogram
from ChromProcess.Processing import integrate_chromatogram_peaks
from ChromProcess.Processing import internal_standard_integral
from ChromProcess.Utils.signal_processing import signal_processing as sig
import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error
chromatogram_directory = f'{experiment_folder}\ChromatogramCSV'
conditions_file = f'{experiment_folder}\{experiment_number}_conditions.csv'
analysis_file = f'{experiment_folder}\{experiment_number}_analysis_details.csv'
peak_collection_directory = f'{experiment_folder}\PeakCollections'


conditions = conditions_from_csv(conditions_file)
analysis = analysis_from_csv(analysis_file)

os.makedirs(peak_collection_directory, exist_ok = True)
chromatogram_files = os.listdir(chromatogram_directory)
chromatogram_files.sort()
chroms = []
for f in chromatogram_files:
    chroms.append(chrom_from_csv(f'{chromatogram_directory}/{f}'))
fit_values = np.array(['amp1','centre1','sigma1','amp2','centre2','sigma2','amp3','centre3','sigma3'])

indices = range(4990,5120)#4820,4871)
fit_error = []
for chrom in chroms:
    chrom_error = []
    lower_bound_test = np.linspace(11.75,11.81,100)
    upper_bound_test = lower_bound_test+0.06
    print('eyo')
    for lb_test, up_test in zip(lower_bound_test,upper_bound_test):
        initial_guess = [15304,lb_test,0.001,15638,11.798,0.00869,17020,11.836,0.009]#[115638,11.45,0.00869]#
        lower_bounds = [0,lb_test,0.0005,10000,11.79,0.005,10000,11.829,0]
        upper_bounds = [1e7,up_test,0.002,50000,11.81,0.012,50000,11.88,0.2]
        #smooth_signal = sig.savitzky_golay(signal, 7, 3, deriv=0, rate=1)
        #chrom.signal=smooth_signal
        popt, pcov = deconvolute_peak(chrom,
                                    initial_guess,
                                    experiment_folder,
                                    indices,
                                    lower_bounds,
                                    upper_bounds,
                                    num_peaks=3,
                                    plotting = False)
        fit_values = np.vstack((fit_values,popt))
        final_curve = np.zeros(len(indices))
        for p in range(0,3):
            gauss = _1gaussian(chrom.time[indices],popt[0+p*3],popt[1+p*3],popt[2+p*3])
            final_curve += gauss
        mse = mean_squared_error(chrom.signal[indices],final_curve)
        chrom_error.append(mse)
        initial_guess=popt
    fit_error.append(chrom_error)
fig, ax = plt.subplots(1,1)
for r in fit_error:
    ax.plot(lower_bound_test,r)
plt.show()

pd.DataFrame(fit_values).to_csv(f'{experiment_folder}\\gaussian_fit_11.70.csv')
pd.DataFrame(fit_error).T.to_csv(f'{experiment_folder}\\error_width.csv')

