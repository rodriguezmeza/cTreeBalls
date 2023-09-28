
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


input1=np.loadtxt("../Output/histZetaM_1.txt")
input2=np.loadtxt("../Output/histZetaM_2.txt")
input3=np.loadtxt("../Output/histZetaM_3.txt")
input4=np.loadtxt("../Output/histZetaM_4.txt")
input5=np.loadtxt("../Output/histZetaM_5.txt")
input6=np.loadtxt("../Output/histZetaM_6.txt")
input7=np.loadtxt("../Output/histZetaM_7.txt")
input8=np.loadtxt("../Output/histZetaM_8.txt")
input9=np.loadtxt("../Output/histZetaM_9.txt")

zetas =[input1/np.mean(input1), input2/np.mean(input2), input3/np.mean(input3),
        input4/np.mean(input4), input5/np.mean(input5), input6/np.mean(input6),
        input7/np.mean(input7), input8/np.mean(input8), input9/np.mean(input9)]
scale = [5, 5, 5, 5, 5, 5, 5, 5, 5]
#fig1 = plt.figure(figsize=(9,6))
#im=plt.imshow(zeta,origin ='lower', cmap='RdBu')
fig1 = plt.figure(figsize=(10, 12))
grid = ImageGrid(fig1, 111, nrows_ncols=(3, 3),
                     axes_pad=0.1,
                     cbar_mode="single",
                     cbar_location="right",)
                     #cbar_size="10%")
for ax, zetasym, scaleplot in zip(grid, zetas, scale):
    im=ax.imshow(zetasym,origin='lower', cmap=cmap, interpolation='bilinear',
                 vmin=-scaleplot, vmax=scaleplot)
#                 norm=MidpointNormalize(midpoint=0), vmin=-scaleplot, vmax=scaleplot)
    #ax.xlabel("$r_1$")
    #plt.ylabel("$r_2$")
#    plt.title("$zeta, m=0$")
#    plt.colorbar(im)

for ax, im_title in zip(grid, ["$m=0$", "$m=1$", "$m=2$", "$m=3$",
                               "$m=4$", "$m=5$", "$m=6$", "$m=7$", "$m=8$", "$m=9$"]):
    t = add_inner_title(ax, im_title, loc='upper left')
    t.patch.set_ec("none")
    t.patch.set_alpha(0.5)

cbar = grid.cbar_axes[0].colorbar(im)
    # Añadir etiquetas a los ejes x e y
for ax in grid:
    ax.set_xlabel("$r_1$")
    ax.set_ylabel("$r_2$")

plt.show()
fig1.savefig(pdf_name_zetaM)

fig5 = plt.figure(figsize=(9,6))
plt.plot(np.diag(input1/np.mean(input1)), label="0")
plt.plot(np.diag(input2/np.mean(input2)), label="1")
plt.plot(np.diag(input3/np.mean(input3)), label="2")
plt.plot(np.diag(input4/np.mean(input4)), label="3")
plt.plot(np.diag(input5/np.mean(input5)), label="4")
plt.plot(np.diag(input6/np.mean(input6)), label="5")
plt.plot(np.diag(input7/np.mean(input7)), label="6")
plt.plot(np.diag(input8/np.mean(input8)), label="7")
plt.plot(np.diag(input9/np.mean(input9)), label="8")
plt.legend()
#plt.xlim(1, 10)
#plt.ylim(-5, 5)
plt.show()
fig5.savefig(pdf_name_diagonal)

