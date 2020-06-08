import ROOT
import array
import math
import os
import numpy as np
import matplotlib.pyplot as plt



def plotHisto(runNumber):
    histolist = [
    "dutdx",
    "sixdxc"    
    ]   
    sigmaDUT = -1000.
    sigmaSIX = -1000.
    runNumber = str(runNumber) 
    fileToOpen = "scopeRD"+runNumber+".root"
    myfile = ROOT.TFile(fileToOpen)

    c1 = ROOT.TCanvas("c1","",900,600)
    c1.SetFillColor(0)
    c1.SetBorderMode(0)
    c1.SetBorderSize(2)
    c1.SetRightMargin(0.15423698)
    c1.SetFrameBorderMode(0)
    c1.SetFrameBorderMode(0) 

    for histoname in histolist:
        myhisto = myfile.Get(histoname)
        try:            
            fitxmin = -0.03
            fitxmax = 0.03
            myhisto.GetXaxis().SetRangeUser(fitxmin, fitxmax)
            myhisto.Draw()
            myFit = myhisto.Fit("gaus","QS","",fitxmin,fitxmax)                
            if "dut" in str(histoname):
                sigmaDUT =myFit.Parameter(2) 
            if "six" in str(histoname):
                sigmaSIX =myFit.Parameter(2) 
        except:
            print histoname, "for run ",runNumber," is missing"
        c1.SaveAs("Run_"+runNumber+"/"+histoname+"_Run"+runNumber+".png")
    if sigmaDUT >0. and sigmaSIX > 0.: 
        deltaRes2 = sigmaDUT*sigmaDUT - sigmaSIX*sigmaSIX/4
        if deltaRes2 > 0:
            print "Run number ",runNumber," sigma DUT ", sigmaDUT, " sigma telescope ",sigmaSIX
            return math.sqrt(deltaRes2)
        else:
            print "ERROR on six/dut resolution!"
            return -1000.
    else:
        return -1000.
#main part
ROOT.gROOT.SetBatch()
ROOT.gErrorIgnoreLevel = ROOT.kFatal
ROOT.gStyle.SetOptStat(0)                                                                                                                                           

runlist = [38826,38827,38828,38829,38830]
angle = [11,10,9,8,6]

sigmax = []
for runNumber in runlist:
    mydir = "Run_"+str(runNumber)
    if not os.path.exists(mydir):
        try:
            os.mkdir(mydir)
        except OSError:
            print ("Creation of the directory %s failed" % path)
    
    sigmax.append(plotHisto(runNumber)*1000.)
                                                                                        

fig =plt.figure(figsize=(7,5), dpi=100)
ax = fig.add_subplot(111)
ax.xaxis.grid(True, which="minor")
ax.yaxis.grid(True, which="major")

plt.xlim(5,14)
#plt.ylim(0,9)
plt.xticks(fontsize = 8)
y_axis_label = "Resolution [um] "
plt.ylabel(y_axis_label)
plt.xlabel("Tilt angle [deg]")
plt.plot(angle,sigmax,'bo',markersize=6,linewidth=0)
ax.patch.set_facecolor("w")
fig.patch.set_facecolor("w")
#adding the legend
from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color="b", lw=3),
                Line2D([0], [0], color="r", lw=3)]

#ax.legend(custom_lines, ['Italy', other_country, 'France'], loc='upper left')
#ax.legend(custom_lines, ['Planar', '3D'], loc='upper right')
plt.show()
