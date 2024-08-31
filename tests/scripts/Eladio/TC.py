import treecorr as tc
import numpy as np
import time
from os import system
from sys import argv


ti = time.time()
print("getting xyz")
tXYZS = time.time()
catfile=np.loadtxt("/pscratch/sd/e/eladio_m/Takahashi_realizations_zs9/cats_zs9_txt/full_sky_whole_XYZK_zs9_r002.txt")

x=catfile.T[0]
y=catfile.T[1]
z=catfile.T[2]
k=catfile.T[3]
cat = tc.Catalog(x=x,y=y,z=z,k=k)
tXYZF = time.time()
print("getting xyz done in {} s".format(np.round(tXYZF-tXYZS,3)))

min_sep = 0.00213811
max_sep = 0.06332045
nbins = 20
max_n = 8

# Lento
#kkk=tc.KKKCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='rad', nbins=nbins,
#                                  max_n=max_n, bin_slop=0, angle_slop=0,
#                                  bin_type='LogMultipole',verbose=2,num_threads=128)

# Rápido
kkk=tc.KKKCorrelation(min_sep=min_sep, max_sep=max_sep, sep_units='rad', nbins=nbins,
                                  max_n=max_n,
                                  bin_type='LogMultipole',verbose=2,num_threads=128)

t1 = time.time()
kkk.process(cat,algo="multipole")  
t2 = time.time()
print('Time for calculating kkk LogMultipole correlation = ',t2-t1)
# Para obtener los multipolos
zeta0=np.real(kkk.zeta[:,:,8]/kkk.ntri[:,:,8])
zeta1=np.real(kkk.zeta[:,:,9]/kkk.ntri[:,:,9])
zeta2=np.real(kkk.zeta[:,:,10]/kkk.ntri[:,:,10])
zeta3=np.real(kkk.zeta[:,:,11]/kkk.ntri[:,:,11])
zeta4=np.real(kkk.zeta[:,:,12]/kkk.ntri[:,:,12])
zeta5=np.real(kkk.zeta[:,:,13]/kkk.ntri[:,:,13])
zeta6=np.real(kkk.zeta[:,:,14]/kkk.ntri[:,:,14])
zeta7=np.real(kkk.zeta[:,:,15]/kkk.ntri[:,:,15])
zeta8=np.real(kkk.zeta[:,:,16]/kkk.ntri[:,:,16])



np.savetxt('zeta_r002.txt',kkk.zeta.flatten())
np.savetxt('ntri_r002.txt',kkk.ntri.flatten())
np.savetxt('rnom1d.txt',kkk.rnom1d)


np.savetxt('zeta0_r002_tc.txt',zeta0)
np.savetxt('zeta1_r002_tc.txt',zeta1)
np.savetxt('zeta2_r002_tc.txt',zeta2)
np.savetxt('zeta3_r002_tc.txt',zeta3)
np.savetxt('zeta4_r002_tc.txt',zeta4)
np.savetxt('zeta5_r002_tc.txt',zeta5)
np.savetxt('zeta6_r002_tc.txt',zeta6)
np.savetxt('zeta7_r002_tc.txt',zeta7)
np.savetxt('zeta8_r002_tc.txt',zeta8)

tf = time.time()
print('The total time is = ',tf-ti)
