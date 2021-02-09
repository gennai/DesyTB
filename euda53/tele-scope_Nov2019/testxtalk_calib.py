import numpy as np 
import matplotlib.pyplot as plt
import sys
runN = str(sys.argv[1])
pippo = np.load("inversedVCAL_w6067_336_793_29_sensormap.npy")
pixelhits = np.genfromtxt("pixelHits_"+runN+".dat",skip_header=1,dtype=int)
#print pippo
q1 = []
q1raw = []
q1rawunpaired = []
q2 = []
q2raw = []
q2rawunpaired = []
#print pixelhits.shape
for i in range(pixelhits.shape[0]):    
    if i%2 ==0 :
        q1.append(pippo[pixelhits[i][0],pixelhits[i][1],pixelhits[i][2]])
        q1raw.append(pixelhits[i][2])        
    else:
        q2.append(pippo[pixelhits[i][0],pixelhits[i][1],pixelhits[i][2]])
        q2raw.append(pixelhits[i][2])

for i in range(pixelhits.shape[0]):    
    if i ==0 : continue
    if i%2 ==0 :
        q1rawunpaired.append(pixelhits[i][2])        
    else:
        q2rawunpaired.append(pixelhits[i][2])


ratioq= []
ratioqraw = []
ratioqrawunpaired = []
for i in range(len(q1)):
    ratioq.append(min(q1[i],q2[i])*1./(q1[i]+q2[i]))
    ratioqraw.append(min(q1raw[i],q2raw[i])*1./(q1raw[i]+q2raw[i]))
for i in range(len(q1rawunpaired)):    
    ratioqrawunpaired.append(min(q1rawunpaired[i],q2rawunpaired[i])*1./(q1rawunpaired[i]+q2rawunpaired[i]))


bins = np.linspace(0, 0.5, 25)
fig, ax = plt.subplots(1,1)
plt.hist(ratioqraw,alpha=0.5,bins=20, label="paired")
plt.hist(ratioqrawunpaired,alpha=0.5,bins=20, label="unpaired")
#plt.hist(ratioq,alpha=1.,bins=50, label="paired calibrated")
plt.legend(loc='upper right')
ax.patch.set_facecolor("w")
fig.patch.set_facecolor("w")
plt.savefig("xtalk_"+runN+".png")
plt.show()
