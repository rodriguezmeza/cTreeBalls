import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import ImageGrid
from IPython import display
import pickle
import matplotlib.pyplot as plt, glob
import scipy as sp


def radtoarc(x):
    return x*180*60/np.pi
######################################
#MODEL
########################################

modelFiles=glob.glob("/pscratch/sd/e/eladio_m/model_zs25/zs25_zetam?.txt") #direccion de los zeta_m para el model
modelFiles.sort()

modelData = []
for i in modelFiles:
    modelData.append(np.loadtxt(i))
modelData=np.array(modelData)

thtArr1=modelData[0].T[0]
thtArr2=modelData[0].T[1]

thetasArr = np.c_[thtArr2,thtArr1]

# for radii, all runs have the same so we only load one of them
rArr_patch = np.loadtxt("/pscratch/sd/e/eladio_m/cTreeBalls/tests/fs_zs25_r000_treeompsincos_binario/rbins.txt") # dirección de los bines

THT1_patch,THT2_patch = np.meshgrid(rArr_patch,rArr_patch)

modelFun=[]
for i in range(len(modelData)):
    modelFun.append(sp.interpolate.NearestNDInterpolator(thetasArr,modelData[i].T[2])(THT1_patch,THT2_patch))
    
    
###################################
# ESTIMADOR
###########################3
n=9 # numero de multipolos

expData0_b = []

y0="/pscratch/sd/e/eladio_m/cTreeBalls/tests/fs_zs25_r000_treeompsincos_binario" #dirección de la carpeta de salida
   
for j in range(1,n+1):
    expData0_b.append(np.loadtxt(y0+"/"+"histZetaM_cos_"+ str(j) +".txt")+np.loadtxt(y0+"/"+"histZetaM_sin_"+ str(j) +".txt"))
    
    
rArr = radtoarc(np.loadtxt(y0+"/rbins.txt")) # bines en arcmin

######################################
# PLOTS
########################################

suptitlefs=24; titlefs=16; lfs=14;tfs = 8; N=len(rArr); stp=4
ls="-";ls2=""; ls3="--"; lw=.9; mrk=""; mrk2="."; mrk3="x";mrk4="+";mrks=6; alph=.6; alph2=.3
for m in [0,1,2,3,4,5]: # para los multipolos
    fig,axs = plt.subplots(nrows=4,ncols=5,figsize=(23,13),sharex=True,sharey=True)
    fig.subplots_adjust(left=0.15,bottom=0.1, right=.85, top=0.9, wspace=0.2, hspace=0.4)
    fig.suptitle(r"$ \zeta_{m =}$"+r"$_{}$".format(m),fontsize=suptitlefs)

    indx = 0 # tht2 index
    nr,nc = np.shape(axs)
    for r in range(nr):
        for c in range(nc):
                # remove frame top and right and set log x- and y-scale
            axs[r,c].spines['top'].set_visible(False) ; axs[r,c].spines['right'].set_visible(False); 
            axs[r,c].set_xscale("log"); axs[r,c].set_yscale("log")

            axs[r,c].plot(rArr,modelFun[m][:,indx],
                            linestyle=ls, marker=mrk, markersize=mrks, color="black", alpha=alph, lw=1,label="model")
            
            axs[r,c].plot(rArr,expData0_b[m][:,indx],
                           linestyle=ls, marker=mrk3, markersize=mrks, color="blue", alpha=alph, lw=1,label="r000")     
                
           
            indx += 1
    fig.text(.7,.95,f"zs25", color="black",fontsize=18)
    handles, labels = axs[nr-1,nc-1].get_legend_handles_labels()
    fig.legend(handles, labels, loc=(.055,.95),frameon=False,fontsize=18,ncols=5)
    #plt.savefig("zs25andmodel_zeta_{}.pdf".format(m))
    
    plt.show()
