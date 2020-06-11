import numpy as np 
pippo = np.load("inversedVCAL_w6067_336_793_29.npy")
ncols = pippo.shape[0]
nrows = pippo.shape[1]
a = np.zeros((ncols/2,nrows*2,16))
col =0
row = 0
for ix in range(ncols):
    for iy in range(nrows):
        col = ix/2
        if ix%2:
	        row = 2*iy + 1
        else:
	        row = 2*iy + 0            
        a[col,row] = pippo[ix,iy]
np.save("inversedVCAL_w6067_336_793_29_sensormap.npy",a)
