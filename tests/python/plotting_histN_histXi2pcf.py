
import numpy as np
from scipy.interpolate import interp2d
from scipy.io import mmread, mmwrite

import matplotlib
import matplotlib.pyplot as plt
#import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import ImageGrid

import mpmath as mp
from mpmath import *

from PIL import Image

#from scipy.spatial import KDTree
#from time import process_time

from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import ImageGrid

cmap=plt.get_cmap('RdBu')
class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

norm = MidpointNormalize(midpoint=0)

def add_inner_title(ax, title, loc, **kwargs):
    from matplotlib.offsetbox import AnchoredText
    from matplotlib.patheffects import withStroke
    prop = dict(path_effects=[withStroke(foreground='w', linewidth=3)],
                size=plt.rcParams['legend.fontsize'])
    at = AnchoredText(title, loc=loc, prop=prop,
                      pad=0., borderpad=0.5,
                      frameon=False, **kwargs)
    ax.add_artist(at)
    return at
    
pdf_name_points = "./points.pdf"
pdf_name_zetaM = "./zetaM.pdf"
pdf_name_diagonal = "./mZeta_diagonal.pdf"
pdf_name_histN = "./histN.pdf"
pdf_name_histXi2pcf = "./histXi2pcf.pdf"


histNFile=np.loadtxt("../Output/histN.txt")
NdataHistN = np.shape(histNFile)[0]
print(NdataHistN)
xd=histNFile[:,0]
yd=histNFile[:,1]
fig6 = plt.figure(figsize=(9,6))
plt.plot(xd,yd)
plt.show()
fig6.savefig(pdf_name_histN)

histNFile=np.loadtxt("../Output/histXi2pcf.txt")
NdataHistN = np.shape(histNFile)[0]
print(NdataHistN)
xd=histNFile[:,0]
yd=histNFile[:,1]
fig7 = plt.figure(figsize=(9,6))
plt.plot(xd,yd)
plt.show()
fig7.savefig(pdf_name_histXi2pcf)


