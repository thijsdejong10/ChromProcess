# from pathlib import Path
import os
import matplotlib.pyplot as plt
from ChromProcess.Loading import conditions_from_csv, chrom_from_csv, chrom_from_cdf
from ChromProcess.Loading.analysis_info.analysis_from_toml import analysis_from_toml
from ChromProcess.Writers import chromatogram_to_csv
from ChromProcess.Processing import ic_background_subtraction
from pathlib import Path

experiment_number = "SPE043"
experiment_folder = Path(f"{Path.home()}/PhD/Data/{experiment_number}")
chromatogram_directory = Path(experiment_folder, f"cdf")
csv_directory = Path(experiment_folder, "MS_chromatogramCSV")
conditions_file = Path(experiment_folder, f"{experiment_number}_conditions.csv")
analysis_file = Path(experiment_folder, f"{experiment_number}_analysis_details.toml")
peak_collection_directory = Path(experiment_folder, f"PeakCollections")
conditions = conditions_from_csv(conditions_file)
analysis = analysis_from_toml(analysis_file)

os.makedirs(csv_directory, exist_ok=True)
chromatogram_files = os.listdir(chromatogram_directory)
chromatogram_files.sort()
chroms = []
for f in chromatogram_files:
    chroms.append(chrom_from_cdf(f"{chromatogram_directory}/{f}"))
for count, c in enumerate(chroms):
    chroms[count].signal = ic_background_subtraction(c, threshold=500)
fig, ax = plt.subplots()
for c in chroms:
    ax.plot(
        c.time[analysis.plot_region[0] : analysis.plot_region[1]],
        c.signal[analysis.plot_region[0] : analysis.plot_region[1]],
        label=c.filename,
    )
plt.show()
for chrom in chroms:
    chromatogram_to_csv(chrom, f"{csv_directory}/{chrom.filename}")
