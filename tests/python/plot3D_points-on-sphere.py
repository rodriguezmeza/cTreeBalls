
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

# Catalog from cute_box:
xyz = np.loadtxt(f"./Output/points_on_sphere.txt")

xs = xyz[:,0]
ys = xyz[:,1]
zs = xyz[:,2]

ax.scatter(xs, ys, zs, marker='o', s=1)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()

fig.savefig(f"data_3D_view.pdf")

