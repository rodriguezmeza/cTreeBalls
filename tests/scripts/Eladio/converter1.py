import numpy as np, healpy as hp, time
from os import system
from sys import argv


fileName = "/pscratch/sd/e/eladio_m/Takahashi_realizations_zs9/Takahashi_sz9/allskymap_nres12r020.zs9.mag.dat"
saveto = "/pscratch/sd/e/eladio_m/Takahashi_realizations_zs9/cats_zs9_txt/"

# for rounding to closest number of decimals, define param C
C = 12

# also use the trasfomation to 3D used in treecorr. ra <-> phi , dec <-> tht 
def getXYZ(ra,dec):
    ra = ra; dec = dec # this expects dec
    if not (len(ra)==len(dec)): return "valError: arrays have diff lenghts"
    # ra = ra; dec = np.pi/2 - dec # this expects collatitude
    lra = len(ra)
    xyz = np.zeros((lra,3))
    
    for i in range(lra):
        # a la treecorr
        xyz[i] = np.cos(dec[i])*np.cos(ra[i]),np.cos(dec[i])*np.sin(ra[i]),np.sin(dec[i])
    return xyz
    
    
print("Start readding from binary")
trs = time.time()
skip = [0, 536870908, 1073741818, 1610612728, 2147483638, 2684354547, 3221225457]       
load_blocks = [skip[i+1]-skip[i] for i in range(0, 6)]
with open(fileName, 'rb') as f:                                                         
    rec = np.fromfile(f, dtype='uint32', count=1)[0]                                    
    nside = np.fromfile(f, dtype='int32', count=1)[0]                                   
    npix = np.fromfile(f, dtype='int64', count=1)[0]                                  
    
    kappa = np.array([])                                                                
    r = npix                                                                            
    for i, l in enumerate(load_blocks):                                                 
        blocks = min(l, r)                                                              
        load = np.fromfile(f, dtype='float32', count=blocks)                            
        np.fromfile(f, dtype='uint32', count=2)                                         
        kappa = np.append(kappa, load)                                                  
        r = r-blocks                                                                    
        if r == 0:                                                                      
            break                                                                       
        elif r > 0 and i == len(load_blocks)-1:                                         
            load = np.fromfile(f, dtype='float32', count=r)                             
            np.fromfile(f, dtype='uint32', count=2)                                     
            kappa = np.append(kappa, load)                                              
        
    # gamma1 = np.array([])                                                               
    # r = npix                                                                            
    # for i, l in enumerate(load_blocks):                                                 
    #     blocks = min(l, r)                                                              
    #     load = np.fromfile(f, dtype='float32', count=blocks)                            
    #     np.fromfile(f, dtype='uint32', count=2)                                         
    #     gamma1 = np.append(gamma1, load)                                                
    #     r = r-blocks                                                                    
    #     if r == 0:                                                                      
    #         break                                                                       
    #     elif r > 0 and i == len(load_blocks)-1:                                         
    #         load = np.fromfile(f, dtype='float32', count=r)                             
    #         np.fromfile(f, dtype='uint32', count=2)                                     
    #         gamma1 = np.append(gamma1, load)                                            
            
    # gamma2 = np.array([])                                                               
    # r = npix                                                                            
    # for i, l in enumerate(load_blocks):                                                 
    #     blocks = min(l, r)                                                              
    #     load = np.fromfile(f, dtype='float32', count=blocks)                            
    #     np.fromfile(f, dtype='uint32', count=2)                                         
    #     gamma2 = np.append(gamma2, load)                                                
    #     r = r-blocks                                                                    
    #     if r == 0:                                                                      
    #         break                                                                       
    #     elif r > 0 and i == len(load_blocks)-1:                                         
    #         load = np.fromfile(f, dtype='float32', count=r) 
    #         np.fromfile(f, dtype='uint32', count=2)                                     
    #         gamma2 = np.append(gamma2, load)                                            
            
    # omega = np.array([])                                                                
    # r = npix                                                                            
    # for i, l in enumerate(load_blocks):                                                 
    #     blocks = min(l, r)                                                              
    #     load = np.fromfile(f, dtype='float32', count=blocks)                            
    #     np.fromfile(f, dtype='uint32', count=2)                                         
    #     omega = np.append(omega, load)                                                  
    #     r = r-blocks                                                                    
    #     if r == 0:                                                                      
    #         break                                                                       
    #     elif r > 0 and i == len(load_blocks)-1:                                         
    #         load = np.fromfile(f, dtype='float32', count=r)                             
    #         np.fromfile(f, dtype='uint32', count=2)                                     
    #         omega = np.append(omega, load)
trf = time.time()            
print('read from binary completed in time {} s\n'.format(np.round(trf-trs,3)))



print("getting angles from kappa array")
tAngDiscS = time.time()
phi, tht = hp.pix2ang(nside = nside, ipix =range(hp.nside2npix(nside)), lonlat = True)
tht = np.deg2rad(tht)
phi = np.deg2rad(phi)
tAngDiscF = time.time()
print("getting angles done in {} s".format(np.round(tAngDiscF-tAngDiscS,3)))



print("getting xyz")
tXYZS = time.time()
xyz = getXYZ(phi, tht)
x = xyz.T[0]
y = xyz.T[1]
z = xyz.T[2]
tXYZF = time.time()
print("getting xyz done in {} s".format(np.round(tXYZF-tXYZS,3)))



print('start to write out realization patches')
tWriteS = time.time()
header1='# nbody NDIM Lx Ly Lz\n'
header2='#'+str(len(x))+' 3 '+str(np.max(x)-np.min(x))+' '+str(np.max(y)-np.min(y))+' '+str(np.max(z)-np.min(z))+'\n'
oFileName = saveto+"full_sky_whole_XYZK_zs9_r020.txt"
with open(oFileName , "w") as f:
    f.write(header1); f.write(header2)
    for j in range(len(x)):
        f.write(str(np.round(x[j],C))+"\t"+str(np.round(y[j],C))\
            +"\t"+str(np.round(z[j],C))+"\t"\
            +str(np.round(kappa[j],C))+"\n")
tWriteF = time.time()
print("write out done in {} s".format(np.round(tWriteF-tWriteS,3)))