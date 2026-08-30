import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(10, 10))

#change path
xi2D = np.loadtxt(f"./Output/histXi2pcf.txt")
xi3D = np.loadtxt(f"./Outputs_to_compare_with/Output_nside256_balls-omp/histXi2pcf.txt")
#

xd1=xi2D[:,0]
yd1=xi2D[:,1]
xd2=xi3D[:,0]
yd2=xi3D[:,1]
fac=180*60/np.pi
ms = 3
idx = 1

plt.plot(xd1*fac, 100*(yd1-yd2)/yd1, linestyle="--", marker="o", label=f"Multipole i={idx}", markersize=ms, c="purple")
plt.xlabel("$\\theta$ [rad]", fontsize=20)
plt.ylabel("$\\theta \\zeta(\\theta)$", fontsize=20);
plt.xscale("log")

# Show plot
#plt.show()

# Save plot
fig.savefig(f"Xi2pcf_err_plots_test_octree-ggg-omp.pdf")

