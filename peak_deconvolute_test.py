import os
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
import numpy as np
import pandas as pd
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
initial_guess = [15304,11.77,0.011066,15638,11.798,0.00869,17020,11.836,0.011]
for chrom in chroms:
    indices = range(4990,5120)
    lower_bounds = [13000,11.76,0.005,10000,11.78,0,10000,11.829,0]
    upper_bounds = [18000,11.78,0.015,30000,11.81,0.2,30000,11.845,0.03]
    #smooth_signal = sig.savitzky_golay(signal, 7, 3, deriv=0, rate=1)
    popt, pcov = deconvolute_peak(chrom,
                                initial_guess,
                                experiment_folder,
                                indices,
                                lower_bounds,
                                upper_bounds,
                                num_peaks=3)
    fit_values = np.vstack((fit_values,popt))

    initial_guess=popt
pd.DataFrame(fit_values).to_csv(f'{experiment_folder}\\gaussian_fit.csv')