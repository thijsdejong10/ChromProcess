import os
from ChromProcess import Classes
from ChromProcess.Loading import peak_collection_from_csv
from ChromProcess.Loading import analysis_from_toml
from ChromProcess.Loading import conditions_from_csv
import numpy as np
from ChromProcess.Loading import chrom_from_csv
experiment_number = 'FRN142'
experiment_folder = r"C:\Users\thijs\Documents\PhD\Data\FRN142"
peak_collection_directory = f'{experiment_folder}\PeakCollections'
conditions_file = f'{experiment_folder}\{experiment_number}_conditions.csv'
analysis_file = f'{experiment_folder}\{experiment_number}_analysis_details.toml'
data_report_directory = f'{experiment_folder}\DataReports'
chromatogram_directory = f'{experiment_folder}\ChromatogramCSV'
os.makedirs(data_report_directory, exist_ok=True)

conditions = conditions_from_csv(conditions_file)
analysis = analysis_from_toml(analysis_file)

chromatogram_files = os.listdir(chromatogram_directory)
chromatogram_files.sort()
chroms = []
for f in chromatogram_files:
    chroms.append(chrom_from_csv(f'{chromatogram_directory}/{f}'))

peak_tables = []
for file in os.listdir(peak_collection_directory):
    if file.endswith('.CSV'):
        peak_tables.append(peak_collection_from_csv(f'{peak_collection_directory}/{file}',round_digits=7))


# Create series of peak collections
series = Classes.PeakCollectionSeries(
                                    peak_tables, 
                                    name = f'{experiment_number}',
                                    conditions = conditions.conditions
                                    )

for pc,chrom in zip(series.peak_collections,chroms):
    for pk in pc.peaks:
        idx = np.argmin(np.abs(pk.retention_time-chrom.time))
        pk.integral= chrom.signal[idx]

IS_pos = 7.08
#series.align_peaks_to_IS(IS_pos)
series.reference_integrals_to_IS()

 # 5% of internal standard integral if integrals are normalised to IS
#series.remove_peaks_below_threshold(peak_removal_limit)
peak_agglomeration_boundary = 0.025 # distance cutoff 

series.get_peak_clusters(bound = peak_agglomeration_boundary)

series.write_data_reports(f'{data_report_directory}/{series.name}', analysis) # create arrays for output
