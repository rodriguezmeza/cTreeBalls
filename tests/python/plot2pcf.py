
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

pdf_name_histXi2pcf = "./histXi2pcf.pdf"

Output_dir = f"./Output/"
extension = ".txt"
histXi2pcfFile_file = f"{Output_dir}histXi2pcf{extension}"

histNFile=np.loadtxt(histXi2pcfFile_file)
NdataHistN = np.shape(histNFile)[0]
print(NdataHistN)
xd=histNFile[:,0]
yd=histNFile[:,1]

fig7 = plt.figure(figsize=(9,6))
plt.plot(xd,yd)
plt.scatter(xd,yd)
plt.xlabel("$\\theta$")
plt.ylabel("$\\zeta$");
plt.show()
fig7.savefig(pdf_name_histXi2pcf)

