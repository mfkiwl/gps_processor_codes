# skyplot.py: Plot class for Fermi observing sky maps
from gbm.data import GbmHealPix
from gbm.plot.gbmplot import SkyHeatmap
import numpy as np
import math
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np

from gbm.data import PosHist
from gbm.plot import SkyPlot
from gbm.plot.gbmplot import SkyPoints # plot element for plotting a single point or many points on the sky

dat = r"/Users/kelemensz/Dropbox/GPS_data/Cumulative_plots/mean_projection/perth/daily_measurements/data_GS_daily_measurements_24h.csv"

def get_matrix(path):
    print('Path: ', path)
    M = pd.read_csv(path, skiprows=0).values#.transpose()[0]#.astype('float64') 
    return M

#============================================================================================================================================================================================================
from gbm.plot.gbmplot import SkyCircle # plot element for plotting a (shaded) circle on the sky

matrix = get_matrix(dat)
ra = np.linspace(-math.pi, math.pi, len(matrix))
dec= np.linspace(-math.pi/2, math.pi/2, len(matrix[0]))
radii = np.random.uniform(1.0, 60.0, size=len(ra)*len(dec))
colors = np.random.random(size=(len(dec),4))

s = SkyPlot( visible=True)
heatmap = SkyHeatmap(ra, dec, matrix.T, s.ax)
heatmap.show()
heatmap.toggle()
plt.show()
# circles = [SkyCircle(ra[i], dec[i], radii[i], s.ax, color=colors[i,:]) for i in range(len(dec))]












