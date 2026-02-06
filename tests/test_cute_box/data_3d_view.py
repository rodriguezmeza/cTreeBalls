#
# Run from tests directory:
#   python test_cute_box/data_3d_view.py
#
import matplotlib.pyplot as plt
import numpy as np


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection='3d')

#plt.style.use('classic')
plt.style.use('bmh')
#plt.style.use('dark_background')
#plt.style.use('fivethirtyeight')
#plt.style.use('ggplot')
#plt.style.use('grayscale')

#xyz = np.loadtxt(f"./Output/points_on_sphere.txt")
#xyz = np.loadtxt(f"./scripts/Abraham/kappa_nres12_zs9NS256r000.txt")
#xyz = np.loadtxt(f"./scripts/Abraham/50patches_zs9_r43/Taka_nres12r043.zs9_region_49of50.txt")
#
# Catalog from cute_box:
xyz = np.loadtxt(f"./test_cute_box/data.txt")
#xyz = np.loadtxt(f"./data.txt")
#xyz = np.loadtxt(f"./test_cute_box/data_format_cballs.txt")

#print(xyz.shape)
xs = xyz[:,0]
ys = xyz[:,1]
zs = xyz[:,2]
#print(xs)

ax.scatter(xs, ys, zs, marker='o', s=1)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()

# Too big (more than 50 Mb)
#fig.savefig(f"data_3D_view.pdf")
