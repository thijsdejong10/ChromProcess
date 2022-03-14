from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt

def plot_peaks_and_boundaries(chroms):
    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})

    for count, chrom in enumerate(chroms):
        ax.plot(chrom.time,[count]*len(chrom.signal),chrom.signal)
        if count%30==0:
            plt.show()
    plt.show()
