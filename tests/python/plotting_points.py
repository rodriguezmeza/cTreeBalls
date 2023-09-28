
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

points=np.loadtxt("../Output/output.txt")
Ndata = np.shape(points)[0]
print(Ndata)
ra=points[:,0]
dec=points[:,1]
kappa=points[:,2]
c = np.sqrt(kappa)
#ra = points[:,1,2]
cat={'ra':ra, 'dec':dec, 'kappa':kappa}
fig1 = plt.figure(figsize=(8, 8))
#plt.scatter(cat['ra'],cat['dec'],marker=".",c=c, s=c, alpha=1.5)
plt.scatter(cat['ra'],cat['dec'], marker=".", s=1)
plt.show()
fig1.savefig(pdf_name_points)


