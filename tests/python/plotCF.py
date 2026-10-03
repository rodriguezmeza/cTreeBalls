
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

pdf_name_histXi2pcf = "./histCF_data.pdf"

normFac=1 #0.85

Output_dir = f"./Output/"
Output_dir_vs = f"./catalogs_working/l-picola/"
#Output_dir_vs = f"./Outputs_to_compare_with/Output_octree-box-omp_halos-0-0-ascii/"
extension = ".txt"
histXi2pcfFile_file = f"{Output_dir}histCF{extension}"
#histXi2pcfFile_file = f"{Output_dir_vs}corr128.dat"
#histXi2pcfFile_file = f"{Output_dir_vs}halos_log.txt"
#histXi2pcfFile_file = f"{Output_dir_vs}halos_nboxes.txt"
histXi2pcfFile_file_vs = f"{Output_dir_vs}halos_nboxes.txt"
#histXi2pcfFile_file_vs = f"{Output_dir_vs}halos_no-I_R_MAX.txt"
#histXi2pcfFile_file_vs = f"{Output_dir_vs}halos_cute-box.txt"
#histXi2pcfFile_file_vs = f"{Output_dir_vs}histCF{extension}"
#histXi2pcfFile_file_vs = f"{Output_dir_vs}corr128{extension}"
#histXi2pcfFile_file_vs = f"{Output_dir_vs}halos{extension}"
#histXi2pcfFile_file_vs = f"{Output_dir_vs}corr128.dat"

histNFile=np.loadtxt(histXi2pcfFile_file)
NdataHistN = np.shape(histNFile)[0]
print(NdataHistN)
xd=histNFile[:,0]
yd=histNFile[:,1]
yd = yd*normFac

histNFile_vs=np.loadtxt(histXi2pcfFile_file_vs)
NdataHistN_vs = np.shape(histNFile_vs)[0]
print(NdataHistN_vs)
xdvs=histNFile_vs[:,0]
ydvs=histNFile_vs[:,1]

fig1 = plt.figure(figsize=(9,6))
#plt.plot(xd,yd)
#plt.scatter(xd,yd)
#plt.plot(xd,xd*yd)
#plt.scatter(xd,xd*yd)
#
plt.plot(xd,xd*xd*yd)
plt.scatter(xd,xd*xd*yd)
#
#plt.plot(xdvs,ydvs)
#plt.plot(xdvs,xdvs*ydvs)
#
plt.plot(xdvs,xdvs*xdvs*ydvs)
#plt.scatter(xdvs,xdvs*xdvs*ydvs)
#
plt.xlabel("$r$")
plt.ylabel("$r \\zeta$");
#plt.xlim(2,100)
#plt.ylim(0.001,10)
#plt.xscale("log")
#plt.yscale("log")
plt.show()
fig1.savefig(pdf_name_histXi2pcf)


fig1 = plt.figure(figsize=(9,6))
plt.plot(xd, 100*(yd-ydvs)/ydvs)
plt.scatter(xd,100*(yd-ydvs)/ydvs)
plt.plot(xd, 100*(ydvs-yd)/yd)
plt.xlabel("$r$")
plt.ylabel("$r \\zeta$");
plt.show()
#fig1.savefig(pdf_name_histXi2pcf)
