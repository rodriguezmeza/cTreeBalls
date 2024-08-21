import numpy as np
import matplotlib.pyplot as plt

# Define la función twopcfplot con ajuste de tamaño de fuente
def zetas(ax, zeta02D, zeta03D, i):
    bins2d = len(zeta2D_0)
    bins3d = len(zeta3D_0)
    zeta02D = [j * (j + 1) * zeta2D_0[j, j + 1] for j in range(1, bins2d-1)]
    zeta03D = [j * (j + 1) * zeta3D_0[j, j + 1] for j in range(1, bins3d-1)]
    ms = 3
    ax.plot(zeta02D, linestyle="--", marker="o", label=f"Region {i} 2D", markersize=ms, c="purple")
    ax.plot(zeta03D, linestyle="--", marker="o", label=f"Region {i} 3D", markersize=ms, c="red")
    ax.set_xlabel(r"$\theta$ [rad]", fontsize=10)
    ax.set_ylabel(r"$j(j+1)\zeta_0$", fontsize=10)
    ax.legend()

# Crear subplots
fig, axs = plt.subplots(10, 5, figsize=(25, 25))

# Iterar a través de los valores de i y cargar los datos
for i, ax in enumerate(axs.flatten()):
    idx = i   # Incrementar índice para saltar valores de i que no se utilizan
    
    if idx == 10 or idx == 14 or idx == 15 or (idx >= 19 and idx <= 20) or (idx >= 24 and idx <= 25) or (idx >= 29 and idx <= 30) or (idx >= 34 and idx <= 35) or (idx >= 39 and idx <= 40):
        continue
    #modificar la ruta
    zeta2D_0 = np.loadtxt(f"./Eladio/Takahashi_tests/takas_043_zs92D/takas2D_{idx}/histZetaM_1.txt")
    zeta3D_0 = np.loadtxt(f"./Eladio/Takahashi_tests/takas_043_zs93D/takas3D_{idx}/histZetaM_1.txt")
    
    zetas(ax, zeta2D_0, zeta3D_0, idx)

# Ajustar el espaciado entre subplots
plt.tight_layout()

# Mostrar los subplots
#plt.show()
fig.savefig(f"zeta_plots.pdf")

