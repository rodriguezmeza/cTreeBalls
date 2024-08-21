#!/usr/bin/env python
# coding: utf-8

import numpy as np
import healpy as hp
from astropy.coordinates import SkyCoord as SC
import glob
import multiprocessing as mp
from os.path import exists
from os import system
import sys

def pprint(nl,totl,token):
    """This function prints each ten percent of a written list
    :nl: number of line
    :totl: total number of lines
    :token: a list containing initially zero (False). Has to be defined outside this function
    """
    
    y = int(np.round((nl/totl)*100, 2))
    if y%10 == 0 and token[0] is False:
        print(f"{nl}: {y} %")
        token[0] = True
    else: token[0] = False

    pass

def write_chunk_to_file(chunk, filename):
    with open(filename, 'a') as f:
        for data_point in chunk:
            line = '\t'.join(map(str, data_point)) + '\n'
            f.write(line)

def write_arrays_to_file_in_parallel(arrays, filename, num_processes=4):
    pool = mp.Pool(processes=num_processes)
    chunk_size = len(arrays[0]) // num_processes

    # Zip arrays together into tuples
    zipped_arrays = list(zip(*arrays))

    # Split zipped arrays into chunks
    chunks = [zipped_arrays[i:i+chunk_size] for i in range(0, len(zipped_arrays), chunk_size)]

    # Write each chunk to the file in parallel
    results = [pool.apply_async(write_chunk_to_file, args=(chunk, filename)) for chunk in chunks]
    pool.close()
    pool.join()

    # Check for any exceptions
    for result in results:
        result.get()


# # Load files from original nside
print("Loading files")

files = glob.glob("./allskymap_nres12r*.zs9.mag.dat_kappa.fits")
files.sort()
print("files to get kappa from",files)

# # Create and export catalog
print("Getting coordinates x y z")

nSide = 4096
nSide2 = nSide
lon, lat = hp.pix2ang(nside = nSide2, ipix = range(hp.nside2npix(nSide2)), lonlat = True)
c = SC(lon, lat, frame="icrs", unit="deg")
XYZ = c.cartesian
x = XYZ.x
y = XYZ.y
z = XYZ.z

print("Export catalogues")

lx = hp.nside2npix(nSide2)
header1 = '# nbody NDIM Lx Ly Lz\n'
header2 = '#' + str(len(x)) + ' 3 ' + str(np.max(x) - np.min(x)) + ' ' + str(np.max(y) - np.min(y)) + ' ' + \
  str(np.max(z) - np.min(z)) + '\n'


def serialExport():

    oPath = f"./T17nres12ToNS{nSide2}_XYZK_prll"
    C = 8
    r = "000 001 002 003 004 005".split()
    names = ["kappa_nres12_zs9r" + ri for ri in r]

    for i,f in enumerate(files):
        d = hp.read_map(f)
        ki = d
        print(f"writting for {names[i]}")
        path = "/global/cfs/cdirs/m1727/2And3PCFCatalogs/wFS/cats/nside4096_zs9_r000005/T17nres12ToNS4096_XYZK_prll"
        with open(path + "/" + names[i] + ".txt", "w") as fo:
            fo.write(header1)
            fo.write(header2)

            printed = [False]
            for j in range(lx):
                fo.write(str(np.round(x[j], C)) + "\t" + str(np.round(y[j],C)) + \
                "\t" + str(np.round(z[j],C)) + "\t" + str(np.round(ki[j],C)) + \
                "\n"
                )
#                pprint(j, lx, printed)

        print(f"done for {names[i]}\n")
    
    pass

def parallelExport():

    C = 8
    #r = "080 081 082 083 084 085 086 087 088 089 090".split()
    r = "000 001 002 003 004 005".split()
    names = [f"kappa_nres12_zs9NS{nSide2}r" + ri for ri in r]

    for i, f in enumerate(files):
        # if i < 7: continue
        # if not exists(f): print(f"file {f} does not exists")
        
        print(f"Loading cat {f}")
        #d = np.loadtxt(f).T[-1] # data.txt
        d = hp.read_map(f)
        
        print("down grading")
        ki = d #hp.ud_grade(d, nSide2)
        
        print(f"writting for {names[i]}")
        oPath = f"./T17nres12ToNS{nSide2}_XYZK_prll"
        system("if ! [[ -f " + oPath + " ]];then mkdir -p " + oPath + " ;fi")
        oFile = oPath + "/" + names[i] + ".txt" 
        with open(oFile, "w") as fo:
            fo.write(header1)
            fo.write(header2)

        data = (x, y, z, ki)
        write_arrays_to_file_in_parallel(data, oFile)

        print(f"Done for {names[i]}")
        
    pass


#parallelExport()
serialExport()
