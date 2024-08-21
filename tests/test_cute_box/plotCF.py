
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

pdf_name_histXi2pcf = "./histCF.pdf"

Output_dir = f"./Output/"
extension = ".txt"
histXi2pcfFile_file = f"{Output_dir}histCF{extension}"

histNFile=np.loadtxt(histXi2pcfFile_file)
NdataHistN = np.shape(histNFile)[0]
print(NdataHistN)
xd=histNFile[:,0]
yd=histNFile[:,1]

fig1 = plt.figure(figsize=(9,6))
plt.plot(xd,yd)
plt.scatter(xd,yd)
plt.xlabel("$r$")
plt.ylabel("$\\zeta$");
plt.show()
fig1.savefig(pdf_name_histXi2pcf)

