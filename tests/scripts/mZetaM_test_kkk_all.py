import numpy as np
import matplotlib.pyplot as plt

maxsubplot = 7
columntoplot = 3
numcolumns = 4
extension = ".txt"
fac=180*60/np.pi
OutputCompare = "./Outputs_to_compare_with/"

#B change path
pathdef1 = "./Output/"
pathdef2 =  OutputCompare + "Output_nside256_balls-omp/"
#pathdef2 =  OutputCompare + "Output_nside256_fits_balls-omp/"
#pathdef2 =  OutputCompare + "Output_nside256_tree-omp-sincos_no-one-ball/"
#pathdef2 =  OutputCompare + "Output_nside256_direct-simple-sincos/"
#pathdef2 =  OutputCompare + "Output_nside256_octree-kkk-omp_NMultipoles_NONORMHIST/"
#pathdef2 =  OutputCompare + "Output_nside1024_balls-omp_zs9r081/"
#pathdef2 =  OutputCompare + "Output_nside2048_balls-omp_zs9r081/"
#pathdef2 =  OutputCompare + "Output_nside2048_octree-kkk-omp_NMultipoles_NONORMHIST_zs9r081/"


search = "comparison"

#E
def inputfile1(type,idx):
    return f"{pathdef1}mhistZetaM_{type}_{idx+1}{extension}"
def inputfile2(type, idx):
    return f"{pathdef2}mhistZetaM_{type}_{idx+1}{extension}"

def savename(type):
    return f"mZetaM_{type}_plots_{search}.pdf"

#B Define plot function with font sizing adjust
def zetas(ax, xmZM1, ymZM1, xmZM2, ymZM2, i):
    ms = 3
    ax.plot(fac*xmZM1, ymZM1, linestyle="--", marker="o", label=f"Multipole i={i}", markersize=ms, c="purple")
    ax.plot(fac*xmZM2, ymZM2, linestyle="None", marker="o", label=f"Multipole i={i} vs", markersize=ms, c="red")
    ax.set_xlabel(r"$\theta$ [rad]", fontsize=10)
    ax.set_ylabel(r"$\zeta_{i}$", fontsize=10)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
#E

#B Define plot function with font sizing adjust
def zetas_err(ax, xmZM1, ymZM1, xmZM2, ymZM2, i):
    ms = 3
    ax.plot(fac*xmZM1, 100*(ymZM1-ymZM2)/ymZM1, linestyle="--", marker="o", label=f"Multipole i={i}", markersize=ms, c="purple")
    ax.set_xlabel(r"$\theta$ [rad]", fontsize=10)
    ax.set_ylabel(r"$\zeta_{i}$", fontsize=10)
    ax.set_xscale("log")
    ax.legend()
#E

# Plotting mZetaM_cos
fig1, axs = plt.subplots(2, numcolumns, figsize=(15, 7.5))
# Loop i and load data
for i, ax in enumerate(axs.flatten()):
    idx = i
    if i > maxsubplot:
        continue
    mZMFile1 = np.loadtxt(inputfile1("cos",idx))
    mZMFile2 = np.loadtxt(inputfile2("cos",idx))
    xmZM1 = mZMFile1[:,0]
    ymZM1 = mZMFile1[:,columntoplot-1]
    xmZM2 = mZMFile2[:,0]
    ymZM2 = mZMFile2[:,columntoplot-1]
    zetas(ax, xmZM1, ymZM1, xmZM2, ymZM2, idx)
# Tight subplots spacing
plt.tight_layout()
# Show subplots
plt.show()
# Save plot
fig1.savefig(savename("cos"))

# Make subplots
fig2, axs = plt.subplots(2, numcolumns, figsize=(15, 7.5))
# Loop i and load data
for i, ax in enumerate(axs.flatten()):
    idx = i
    if i > maxsubplot:
        continue
    mZMFile1 = np.loadtxt(inputfile1("sin",idx))
    mZMFile2 = np.loadtxt(inputfile2("sin",idx))
    xmZM1 = mZMFile1[:,0]
    ymZM1 = mZMFile1[:,columntoplot-1]
    xmZM2 = mZMFile2[:,0]
    ymZM2 = mZMFile2[:,columntoplot-1]
    zetas(ax, xmZM1, ymZM1, xmZM2, ymZM2, idx)
# Tight subplots spacing
plt.tight_layout()
# Show subplots
plt.show()
# Save plot
fig2.savefig(savename("sin"))

# Make subplots
fig3, axs = plt.subplots(2, numcolumns, figsize=(15, 7.5))
for i, ax in enumerate(axs.flatten()):
    idx = i
    if i > maxsubplot:
        continue
    mZMFile1 = np.loadtxt(inputfile1("cos",idx))
    mZMFile2 = np.loadtxt(inputfile2("cos",idx))
# Loop i and load data
for i, ax in enumerate(axs.flatten()):
    idx = i
    if i > maxsubplot:
        continue
    mZMFile1cos = np.loadtxt(inputfile1("cos",idx))
    mZMFile1sin = np.loadtxt(inputfile1("sin",idx))
    mZMFile2cos = np.loadtxt(inputfile2("cos",idx))
    mZMFile2sin = np.loadtxt(inputfile2("sin",idx))
    xmZM1 = mZMFile1[:,0]
    ymZM1 = mZMFile1cos[:,columntoplot-1] + mZMFile1sin[:,columntoplot-1]
#    ymZM1 = mZMFile1cos[:,columntoplot-1]
    xmZM2 = mZMFile2[:,0]
    ymZM2 = mZMFile2cos[:,columntoplot-1] + mZMFile2sin[:,columntoplot-1]
    zetas(ax, xmZM1, ymZM1, xmZM2, ymZM2, idx)
# Tight subplots spacing
plt.tight_layout()
# Show subplots
plt.show()
# Save plot
fig3.savefig(savename("cos+sin"))

# Make err plots
# Make subplots
fig4, axs = plt.subplots(2, numcolumns, figsize=(15, 7.5))
for i, ax in enumerate(axs.flatten()):
    idx = i
    if i > maxsubplot:
        continue
    mZMFile1 = np.loadtxt(inputfile1("cos",idx))
    mZMFile2 = np.loadtxt(inputfile2("cos",idx))
# Loop i and load data
for i, ax in enumerate(axs.flatten()):
    idx = i
    if i > maxsubplot:
        continue
    mZMFile1cos = np.loadtxt(inputfile1("cos",idx))
    mZMFile1sin = np.loadtxt(inputfile1("sin",idx))
    mZMFile2cos = np.loadtxt(inputfile2("cos",idx))
    mZMFile2sin = np.loadtxt(inputfile2("sin",idx))
    xmZM1 = mZMFile1[:,0]
    ymZM1 = mZMFile1cos[:,columntoplot-1] + mZMFile1sin[:,columntoplot-1]
#    ymZM1 = mZMFile1cos[:,columntoplot-1]
    xmZM2 = mZMFile2[:,0]
    ymZM2 = mZMFile2cos[:,columntoplot-1] + mZMFile2sin[:,columntoplot-1]
    zetas_err(ax, xmZM1, ymZM1, xmZM2, ymZM2, idx)
# Tight subplots spacing
plt.tight_layout()
# Show subplots
plt.show()
# Save plot
fig4.savefig(savename("cos+sin_err"))
