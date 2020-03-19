import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

deltaRays = pd.read_csv('./deltaRays.csv', index_col=0)
deltaRays = pd.DataFrame(deltaRays)

cols = deltaRays.loc[:," col"]
rows = deltaRays.loc[:," row"]
indeces = deltaRays.index
subcols = []
subrows = []
subcolors= []
subindices = []
j=0
for i in range(20000):
    subcols.append(cols.iloc[i])
    subrows.append(rows.iloc[i])
    subindices.append(indeces[i])
    if(indeces[i] == 0):
        j=j+1
        color = j*2./3000
    subcolors.append(color)
    #print(cols.iloc[i], rows.iloc[i], indeces[i])
#colors = cm.rainbow(np.linspace(0, 1, len(subcols)))
colors = cm.rainbow(subcolors)
#for a,i,j,k in zip(subindices,subcols,subrows,subcolors):
#    print(a,i,j,k)
ax = plt.figure(figsize=(6,8), dpi=100).add_subplot(111)
plt.scatter(subcols,subrows, color=colors, s=2.0)
plt.xlabel("Columns")
plt.ylabel("Rows")
plt.show()