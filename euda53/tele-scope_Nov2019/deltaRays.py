import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

deltaRays = pd.read_csv('./deltaRays.csv', index_col=0)
deltaRays = pd.DataFrame(deltaRays)

cols = deltaRays.loc[:,"col"]
xpos  = deltaRays.loc[:,"xpos"]
ypos  = deltaRays.loc[:,"ypos"]
rows = deltaRays.loc[:,"row"]
tots = deltaRays.loc[:,"tot"]
indeces = deltaRays.index
subcols = []
subrows = []
subcolors= []
subindices = []
subxpos = []
subypos = []
j=0
size = []

#print(max(rows))
for i in range(2000):
    subcols.append(cols.iloc[i])
    subrows.append(rows.iloc[i])
    subindices.append(indeces[i])
    subxpos.append(xpos.iloc[i])
    subypos.append(ypos.iloc[i])
    size.append(1.*tots.iloc[i])
    #print(indeces[i], tots.iloc[i], cols.iloc[i],rows.iloc[i] )
    if(indeces[i] == 0):
        j=j+1
        color = j*3./200
    subcolors.append(color)
    #print(cols.iloc[i], rows.iloc[i], indeces[i])
#colors = cm.rainbow(np.linspace(0, 1, len(subcols)))
colors = cm.rainbow(subcolors)
print (colors)
#for a,i,j,k in zip(subindices,subcols,subrows,subcolors):
#    print(a,i,j,k)
ax = plt.figure(figsize=(6,8), dpi=100).add_subplot(111)
plt.scatter(subcols,subrows, color=colors, s=size)
plt.scatter(subxpos,subypos, color=colors,marker="x")
plt.xlabel("Columns")
plt.ylabel("Rows")
plt.show()