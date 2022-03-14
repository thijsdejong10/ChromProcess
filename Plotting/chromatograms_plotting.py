from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
def plot_peaks_and_boundaries(chroms):
    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
    z = []
    X,Y = np.meshgrid(chroms[0].time[:12990],np.linspace(0,len(chroms),len(chroms)))
    for count, chrom in enumerate(chroms):
        z.append(chrom.signal[20:12990])
    
    z = np.array(z)
    ax.plot_wireframe(X,Y,z)
    plt.show()

def heatmap_cluster(chroms):
    from ChromProcess.Utils.utils.clustering import cluster
    peak_pos = np.array([])
    for chrom in chroms:
        peak_pos = np.hstack((peak_pos, np.array(chrom.peaks.keys())))
    peaks = np.argsort(peak_pos)
    clusts = []
    for c in cluster(peaks, bound=0.025):
        clusts.append(c)

    z = []
    for count, chrom in enumerate(chroms):
        z.append(chrom.signal[20:12990])
    
    z = np.array(z)
    fig = plt.figure()
    extent = [chroms[0].time[0], chroms[0].time[12990], 0, len(chroms)]
    plt.imshow(np.log(z+0.001),extent=extent,aspect=0.01)
    


    plt.show()
