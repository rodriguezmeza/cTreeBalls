
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

pdf_name_histXi2pcf = "./histCF_cola_064.pdf"

Output_dir = f"./Output/"
Output_dir_vs = f"./test_cute_box/output_cute_box/"
extension = ".txt"
extension_vs = ".dat"
histXi2pcfFile_file = f"{Output_dir}histCF{extension}"
histXi2pcfFile_file_vs = f"{Output_dir_vs}corr128_cola_064{extension_vs}"

histNFile=np.loadtxt(histXi2pcfFile_file)
NdataHistN = np.shape(histNFile)[0]
print(NdataHistN)
xd=histNFile[:,0]
yd=histNFile[:,1]

histNFile_vs=np.loadtxt(histXi2pcfFile_file_vs)
NdataHistN_vs = np.shape(histNFile_vs)[0]
print(NdataHistN_vs)
xdvs=histNFile_vs[:,0]
ydvs=histNFile_vs[:,1]

fig1 = plt.figure(figsize=(9,6))
plt.plot(xd,xd*yd)
plt.scatter(xd,xd*yd)
plt.plot(xdvs,xdvs*ydvs)
plt.xlabel("$r$")
plt.ylabel("$r \\zeta$");
plt.show()
fig1.savefig(pdf_name_histXi2pcf)

