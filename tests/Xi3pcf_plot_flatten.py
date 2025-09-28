import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec

n=8 # número de multipolos

#B changes here
y1 = "./Output/"

y1fileCos = y1+"/"+"histZetaM_cos_"
y1fileSin = y1+"/"+"histZetaM_sin_"
ofilename="flatten.pdf"
y1Xi2pcf=y1+"histXi2pcf.txt"
ofilenameXi2pcf=f"Xi2pcf.pdf"

xlim=400
ylim1=-1e-8
ylim2=1e-8
#E

xlabel = r"triangle index"
ylabel0 = r"$ \zeta_0$"
ylabel1 = r"$ \zeta_1$"
ylabel2 = r"$ \zeta_2$"
ylabel3 = r"$ \zeta_3$"
ylabel4 = r"$ \zeta_4$"
xlabelXi2pcf = "$\\theta$ [rad]"
ylabelXi2pcf = "$\\theta \\zeta(\\theta)$"

label1 = "case 1"

color0 = "blue"

symbolsize = 5
lw = 1.75
#
linestyle1 = "-"
#

expData1 = []

for j in range(1,n+1):
    expData1.append(np.loadtxt(y1fileCos + str(j) +".txt") + np.loadtxt(y1fileSin + str(j) +".txt"))
#

fig = plt.figure(tight_layout=True,figsize=(15, 15))
gs = gridspec.GridSpec(5, 3)
formatter=mpl.ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)

# m=0
ax = fig.add_subplot(gs[0, :])
ax.plot(range(len(expData1[0].flatten())),expData1[0].flatten(),color=color0,linestyle=linestyle1,linewidth=lw,alpha=1,label=label1)
ax.set_ylabel(ylabel0,fontsize=35)
ax.set_xlim(0,xlim)
ax.tick_params(labelsize = 25)
#ax.legend(fontsize=25)
ax.axes.xaxis.set_ticklabels([])
# m=1
ax = fig.add_subplot(gs[1, :])
ax.plot(range(len(expData1[1].flatten())),expData1[1].flatten(),color=color0,linestyle=linestyle1,linewidth=lw,alpha=1)
ax.set_ylabel(ylabel1,fontsize=35)
ax.set_xlim(0,xlim)
ax.tick_params(labelsize = 25)
#ax.legend(fontsize=25)
ax.axes.xaxis.set_ticklabels([])
# m=2
ax = fig.add_subplot(gs[2, :])
ax.plot(range(len(expData1[2].flatten())),expData1[2].flatten(),color=color0,linestyle=linestyle1,linewidth=lw,alpha=1)
ax.set_ylabel(ylabel2,fontsize=35)
ax.set_xlim(0,xlim)
ax.tick_params(labelsize = 25)
#ax.legend(fontsize=25)
# m=3
ax = fig.add_subplot(gs[3, :])
ax.plot(range(len(expData1[3].flatten())),expData1[3].flatten(),color=color0,linestyle=linestyle1,linewidth=lw,alpha=1)
ax.set_ylabel(ylabel3,fontsize=35)
ax.set_xlim(0,xlim)
ax.tick_params(labelsize = 25)
#ax.legend(fontsize=25)
# m=4
ax = fig.add_subplot(gs[4, :])
ax.plot(range(len(expData1[4].flatten())),expData1[4].flatten(),color=color0,linestyle=linestyle1,linewidth=lw,alpha=1)
ax.set_xlabel(xlabel,fontsize=30)
ax.set_ylabel(ylabel4,fontsize=35)
ax.set_xlim(0,xlim)
ax.tick_params(labelsize = 25)
#ax.legend(fontsize=25)

fig.align_labels()
plt.savefig(ofilename,dpi=300)
#plt.show()


#B Plotting 2pcf Xi
xi1D = np.loadtxt(y1Xi2pcf)
xd1=xi1D[:,0]
yd1=xi1D[:,1]
fac=180*60/np.pi
ms = 3
fig = plt.figure(figsize=(10, 10))
plt.plot(xd1*fac,fac*xd1*yd1, linestyle="--", marker="o", label=label1, markersize=ms, c="purple")
plt.xlabel(xlabelXi2pcf, fontsize=20)
plt.ylabel(ylabelXi2pcf, fontsize=20);
# Show plot
#plt.show()
# Save plot
fig.savefig(ofilenameXi2pcf)
#E


