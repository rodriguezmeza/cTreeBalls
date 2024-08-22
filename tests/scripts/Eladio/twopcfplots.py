import numpy as np
import matplotlib.pyplot as plt


def twopcfplot(ax, xi2D, xi3D, i):
    ms = 3
    ax.plot(xi2D[:, 0], xi2D[:, 0] * xi2D[:, 1]*2, linestyle="--", marker="o", label=f"Region {i} 2D", markersize=ms, c="purple")
    ax.plot(xi3D[:, 0], xi3D[:, 0] * xi3D[:, 1]*2, linestyle="--", marker="o", label=f"Region {i} 3D", markersize=ms, c="red")
    ax.set_xlabel(r"$\theta$ [rad]", fontsize=10)
    ax.set_ylabel(r"$\theta \, \xi_{kk}$", fontsize=10)
    ax.legend()


fig, axs = plt.subplots(10, 5, figsize=(30, 30))


for i, ax in enumerate(axs.flatten()):
    idx = i   
    
    if idx == 10 or idx == 14 or idx == 15 or (idx >= 19 and idx <= 20) or (idx >= 24 and idx <= 25) or (idx >= 29 and idx <= 30) or (idx >= 34 and idx <= 35) or (idx >= 39 and idx <= 40):
        continue
    #modificar la ruta
    xi2D = np.loadtxt(f"./Eladio/Takahashi_tests/takas_043_zs92D/takas2D_{idx}/histXi2pcf.txt")
    xi3D = np.loadtxt(f"./Eladio/Takahashi_tests/takas_043_zs93D/takas3D_{idx}/histXi2pcf.txt")
    
    twopcfplot(ax, xi2D, xi3D, idx)


plt.tight_layout()


#plt.show()
fig.savefig(f"twopcf_plot.pdf")

