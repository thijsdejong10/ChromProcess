#from pathlib import Path
import os
import matplotlib.pyplot as plt
from ChromProcess.Loading import chrom_from_csv, experiment_conditions
from ChromProcess.Loading import analysis_from_csv
from ChromProcess.Loading import conditions_from_csv
from ChromProcess.Loading.analysis_info.analysis_from_toml import analysis_from_toml
from ChromProcess.Loading.peak.peak_from_csv import peak_indices_from_file, peak_boundaries_from_file
from ChromProcess.Utils.signal_processing.deconvolution import deconvolute_peak
experiment_number = 'FRN146'
experiment_folder = f"C:/Users/thijs/Documents/PhD/Data/{experiment_number}"
from ChromProcess.Utils.peak_finding import find_peaks_scipy
from ChromProcess.Utils import indices_from_boundary, peak_indices_to_times
import pandas as pd
from ChromProcess.Processing import add_peaks_to_chromatogram
from ChromProcess.Processing import integrate_chromatogram_peaks
from ChromProcess.Processing import internal_standard_integral_look_ahead
import numpy as np
from ChromProcess import Classes
from Plotting.chromatograms_plotting import plot_peaks_and_boundaries, heatmap_cluster, peak_area

chromatogram_directory = f'{experiment_folder}/ChromatogramCSV'
conditions_file = f'{experiment_folder}/{experiment_number}_conditions.csv'
analysis_file = f'{experiment_folder}/{experiment_number}_analysis_details.toml'
peak_collection_directory = f'{experiment_folder}/PeakCollections'

conditions = conditions_from_csv(conditions_file)
analysis = analysis_from_toml(analysis_file)

os.makedirs(peak_collection_directory, exist_ok = True)
chromatogram_files = os.listdir(chromatogram_directory)
chromatogram_files.sort()
chroms = []
for f in chromatogram_files:
    chroms.append(chrom_from_csv(f'{chromatogram_directory}/{f}'))
#blank = chrom_from_csv(f'{experiment_folder}\\control_injections\\FRN142_blank_2.dx_FID1A.CSV')
#subtract baseline
#for chrom in chroms:
#    chrom.signal -= blank.signal


# fig, ax = plt.subplots()
# for c in chroms:
#     ax.plot(c.time[analysis.plot_region[0]:analysis.plot_region[1]], 
#             c.signal[analysis.plot_region[0]:analysis.plot_region[1]],
#             label = c.filename)
# plt.show()
is_start = analysis.internal_standard_region[0]
is_end = analysis.internal_standard_region[1]
for c in chroms:
    c.signal = c.signal - min(c.signal[analysis.plot_region[0]:analysis.plot_region[1]])
    internal_standard_integral_look_ahead(c, is_start, is_end)
    c.signal = c.signal/c.internal_standard.integral
    # internal_standard_integral_look_ahead(c, is_start, is_end)
    # print(c.internal_standard.integral)

fig, ax = plt.subplots()
for c in chroms:
    ax.plot(c.time[analysis.plot_region[0]:analysis.plot_region[1]], 
            c.signal[analysis.plot_region[0]:analysis.plot_region[1]],
            label = c.filename)
plt.show()

threshold = analysis.peak_pick_threshold
if type(threshold) == float:
    threshold = [threshold for r in analysis.regions]
peak_figure_folder = f"{experiment_folder}\\peak_figures"
os.makedirs(peak_figure_folder,exist_ok=True)
for chrom in chroms:
    for reg,thres in zip(analysis.regions,threshold):
        inds = indices_from_boundary(chrom.time, reg[0], reg[1])
        time = chrom.time[inds]
        signal = chrom.signal[inds]
        picked_peaks = find_peaks_scipy(signal, 
                        threshold=thres, 
                        min_dist=analysis.peak_distance, 
                        max_inten = 1e100, 
                        prominence = analysis.prominence, 
                        wlen = 1001, 
                        look_ahead = analysis.boundary_window,
                        smooth_window=41)
        peak_features = peak_indices_to_times(time,picked_peaks)
        peaks = []
        for x in range(0, len(picked_peaks["Peak_indices"])):
            pk_idx = picked_peaks["Peak_indices"][x]
            start_idx = picked_peaks["Peak_start_indices"][x]
            end_idx = picked_peaks["Peak_end_indices"][x]
    
            retention_time = time[pk_idx]
            start = time[start_idx]
            end = time[end_idx]
            height = signal[pk_idx] - min(signal) #subtract the baseline of the region from the peak height
            peaks.append(Classes.Peak(retention_time, start, end, indices=[], height=height ))
        peak_area(time,signal,picked_peaks,save_folder=f"{peak_figure_folder}\\{reg[0]}_{chrom.filename[:-4]}.png")
        add_peaks_to_chromatogram(peaks, chrom)
    integrate_chromatogram_peaks(chrom,baseline_subtract=True)


#heatmap_cluster(chroms)
for reg in analysis.deconvolve_regions:
    region_start = analysis.deconvolve_regions[reg]["region_boundaries"][0]
    indices = indices_from_boundary(
            chrom.time, 
            analysis.deconvolve_regions[reg]["region_boundaries"][0],
            analysis.deconvolve_regions[reg]["region_boundaries"][1]
            )
    peak_folder = f"{experiment_folder}\\deconvolved_peaks\\{region_start}"
    os.makedirs(peak_folder,exist_ok=True)
    fit_values = np.array(['mse'])
    for n in range(1,analysis.deconvolve_regions[reg]["number_of_peaks"]+1):
        fit_values = np.hstack((fit_values, [f'amp{n}',f'centre{n}',f'sigma{n}']))
    fit_values = np.hstack((fit_values,'baseline'))
    for chrom in chroms:
        try: #the fit can often crash because of boundaries that do not match, this is costly as it means all the code has to be run again and we don't convolute all regions, try except here makes this easier.
            popt, pcov, mse, peaks = deconvolute_peak(chrom,
                                peak_folder,
                                indices,
                                analysis.deconvolve_regions[reg],
                                plotting = True)
        except:
            print(f"error in fitting {reg}")
        fit_values = np.vstack((fit_values,np.array([mse,*popt])))
        k = [*chrom.peaks.keys()]
        v = [*chrom.peaks.values()]
        for peak in peaks:
            rt = peak.retention_time
            idx = np.where(
                (chrom.time >= peak.start) & (chrom.time <= peak.end)
            )[0]
            peak.indices = idx
            insert = np.searchsorted(k,rt)
            k.insert(insert,rt)
            v.insert(insert,peak)
        chrom.peaks = dict(zip(k,v))   
        
    pd.DataFrame(fit_values).to_csv(f'{peak_folder}\\gaussian_fit_{region_start}.csv')
#for chrom in chroms:
#    peaks_indices = peak_indices_from_file(chrom,f"{peak_collection_directory}\\{chrom.filename}")
#    peak_starts, peak_ends = peak_boundaries_from_file(chrom,f"{peak_collection_directory}\\{chrom.filename}")
#    picked_peaks = {'Peak_indices':peaks_indices, 'Peak_start_indices':peak_starts, 'Peak_end_indices':peak_ends}
#    
#    peak_features = peak_indices_to_times(chrom.time,picked_peaks)
#    add_peaks_to_chromatogram(peak_features, chrom)
#    integrate_chromatogram_peaks(chrom)

#print('test')


heatmap_cluster(chroms,analysis.plot_region)
for c,v in zip(chroms, conditions.series_values):
    c.write_peak_collection(filename = f'{peak_collection_directory}/{c.filename}',
                        header_text = f"{conditions.series_unit},{v}\n",
                        )
