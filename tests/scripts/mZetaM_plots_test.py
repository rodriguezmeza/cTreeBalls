import numpy as np
import matplotlib.pyplot as plt

maxsubplot = 7
columntoplot = 3
extension = ".txt"
fac=180*60/np.pi

#B change path
pathdef1 = "./Output/"
#pathdef2 = "./Outputs_to_compare_with/Output_nside256_balls-omp/"
#pathdef2 = "./Outputs_to_compare_with/Output_nside256_tree-omp-sincos/"
#pathdef2 = "./Outputs_to_compare_with/Output_nside256_tree-omp-sincos_no-one-ball/"
pathdef2 = "./Outputs_to_compare_with/Output_nside256_direct-sincos/"
#E
def inputfile1(type,idx):
    return f"{pathdef1}mhistZetaM_{type}_{idx+1}{extension}"
def inputfile2(type, idx):
    return f"{pathdef2}mhistZetaM_{type}_{idx+1}{extension}"

def savename(type):
    return f"mZetaM_{type}_plots.pdf"

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


# Make plots

# Plotting histNN
histNFile1 = np.loadtxt(f"{pathdef1}histNN{extension}")
histNFile2 = np.loadtxt(f"{pathdef2}histNN{extension}")
xd1=histNFile1[:,0]
yd1=histNFile1[:,1]
xd2=histNFile2[:,0]
yd2=histNFile2[:,1]
fig0 = plt.figure(figsize=(9,6))
plt.plot(fac*xd1,yd1)
plt.scatter(fac*xd1,yd1)
plt.plot(fac*xd2,yd2)
plt.scatter(fac*xd2,yd2)
plt.show()
fig0.savefig(savename("histNN"))

# Plotting CF
histNFile1 = np.loadtxt(f"{pathdef1}histCF{extension}")
histNFile2 = np.loadtxt(f"{pathdef2}histCF{extension}")
xd1=histNFile1[:,0]
yd1=histNFile1[:,1]
xd2=histNFile2[:,0]
yd2=histNFile2[:,1]
fig0 = plt.figure(figsize=(9,6))
plt.plot(fac*xd1,yd1)
plt.scatter(fac*xd1,yd1)
plt.plot(fac*xd2,yd2)
plt.scatter(fac*xd2,yd2)
plt.show()
fig0.savefig(savename("histCF"))

# Plotting Xi2pcf
histNFile1 = np.loadtxt(f"{pathdef1}histXi2pcf{extension}")
histNFile2 = np.loadtxt(f"{pathdef2}histXi2pcf{extension}")
xd1=histNFile1[:,0]
yd1=histNFile1[:,1]
xd2=histNFile2[:,0]
yd2=histNFile2[:,1]
fig0 = plt.figure(figsize=(9,6))
plt.plot(fac*xd1,yd1)
plt.scatter(fac*xd1,yd1)
plt.plot(fac*xd2,yd2)
plt.scatter(fac*xd2,yd2)
plt.show()
fig0.savefig(savename("histXi2pcf"))

fig0 = plt.figure(figsize=(9,6))
plt.plot(fac*xd1,100*(yd1-yd2)/yd1)
plt.scatter(fac*xd1,100*(yd1-yd2)/yd1)
plt.xscale("log")
plt.show()
fig0.savefig(savename("histXi2pcf_err"))

# Plotting mZetaM_cos
fig1, axs = plt.subplots(2, 5, figsize=(15, 7.5))
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
fig2, axs = plt.subplots(2, 5, figsize=(15, 7.5))
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
fig3, axs = plt.subplots(2, 5, figsize=(15, 7.5))
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
    xmZM2 = mZMFile2[:,0]
    ymZM2 = mZMFile2cos[:,columntoplot-1] + mZMFile2sin[:,columntoplot-1]
    zetas(ax, xmZM1, ymZM1, xmZM2, ymZM2, idx)
# Tight subplots spacing
plt.tight_layout()
# Show subplots
plt.show()
# Save plot
fig3.savefig(savename("cos+sin"))

# Make subplots
fig4, axs = plt.subplots(2, 5, figsize=(15, 7.5))
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
    xmZM2 = mZMFile2[:,0]
    ymZM2 = mZMFile2cos[:,columntoplot-1] + mZMFile2sin[:,columntoplot-1]
    zetas_err(ax, xmZM1, ymZM1, xmZM2, ymZM2, idx)
# Tight subplots spacing
plt.tight_layout()
# Show subplots
plt.show()
# Save plot
fig4.savefig(savename("cos+sin_err"))
