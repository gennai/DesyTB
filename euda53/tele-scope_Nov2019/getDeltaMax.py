import ROOT
import numpy as np
import matplotlib.pyplot as plt
import sys

histoname = str(sys.argv[1])
runlist = [39313,39330,39324,39335,39334,39333,39321]
vThr = [1150,1050,1050,1400,1150,950,1315]
vbias = [200,800,500,500,500,500,500]
ratios = []
for (thr, run, bias) in zip(vThr,runlist, vbias):
    fileToOpen = "scopeRD"
    fileToOpen += str(run)+".root"
    print fileToOpen
    myfile = ROOT.TFile(fileToOpen)
    myhistobond =   myfile.Get(histoname+"bonding")
    myhistonobond = myfile.Get(histoname+"nobonding")
    ratio = myhistobond.GetMaximum() - myhistobond.GetMinimum()
    ratio = ratio / (myhistonobond.GetMaximum() - myhistonobond.GetMinimum())    
    ratios.append(ratio)
    print " Run ",run," Vthr ", thr, " Vbias ",bias," ratio ", ratio, \
          " delta bond ", round(myhistobond.GetMaximum() - myhistobond.GetMinimum(),4), \
          " delta no bond ", round(myhistonobond.GetMaximum() - myhistonobond.GetMinimum(),4)

fig =plt.figure(figsize=(7,5), dpi=100)
ax = fig.add_subplot(111)
ax.xaxis.grid(True, which="minor")
ax.yaxis.grid(True, which="major")
plt.xlim(900,1500)
plt.ylim(1.5,3.)
plt.xticks(fontsize = 8)
y_axis_label = " Ratio = Delta Bond / Delta NoBond"
plt.ylabel(y_axis_label)
plt.xlabel("Threshold [electrons]")
plt.plot(vThr,ratios,'bo',markersize=6,linewidth=0, label="Vbias 500 V")
plt.plot(vThr[0],ratios[0],'go',markersize=6,linewidth=0, label = "Vbias 200 V")
plt.plot(vThr[1],ratios[1],'ro',markersize=6,linewidth=0, label = "Vbias 800 V")
plt.legend(loc='upper left', numpoints =1)
#plt.plot(vthr[0],ratios[0],'ro',markersize=6,linewidth=0)
ax.patch.set_facecolor("w")
fig.patch.set_facecolor("w")
#adding the legend
#plt.show()
plt.savefig("ratioVsThreshold_"+histoname+".pdf")
