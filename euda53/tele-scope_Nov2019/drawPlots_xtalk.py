import ROOT
import array
import math
import os
import numpy as np
import matplotlib.pyplot as plt


def plotHisto(runNumber):
    histolist = [
    "linq",
    "effvsxmym",
    "effvsxy",
#    "dutpxqvsxy",
    "linnpxvsxmym",
    "linqxvsxmym",
#    "linqxvsxmymaverage",
    "linqxvsxy",
    "linq",
    "linrowmin2",
    "dutdxdy",
#    "dutdyc2",
    "linqMinOverqCluster",
    "effvst3",
    "linnrow1odd",
    "sixdtx",
    "sixdty",
    "sixdtyLargeClusters", #slope on x for driplets
    "sixdtxLargeClusters", #slope on y for driplets
    "linnrow1eve"
    ]   
    
    x = [-3.5,-3.5,3.0,3.0,-3.5]
    y = [-4.6,4.6,4.6,-4.6,-4.6]


    xp = xbins = array.array('d',x)
    yp = xbins = array.array('d',y)
    

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
    oddEntries = 0.
    evenEntries = 0.
    ratioNumbers = []
    meanFromMoyalHisto = 0 
    meanFromMoyalHistoWithMoyalWidth = 0         
    for histoname in histolist:
        myhisto = myfile.Get(histoname)
        try:
            #if ("sixdt" in histoname):
            #    print histoname,"\n"
            #    print "mean and rms ", round(myhisto.GetMean(),6), round(myhisto.GetRMS(),6),"\n"
            if(histoname =="linnrow1odd"):
                oddEntries = myhisto.GetBinContent(2)
            if(histoname =="linnrow1eve"):
                evenEntries = myhisto.GetBinContent(2)    
            if (histoname =="linqMinOverqCluster"):
                myhisto.Rebin(2)
                myhisto.SetMaximum(35000)
            if (histoname =="linrowmin2"):
                myhisto.GetXaxis().SetRangeUser(100.,140.)
                myhisto.Draw()
                c1.SaveAs("Run_"+runNumber+"/"+histoname+"_Run"+runNumber+".pdf")              
            if (histoname != "effvst3" and histoname != "linq" and histoname != "linrowmin2"):                                                  
                myhisto.Draw("colz")
                myhisto.GetZaxis().SetTitleOffset(1.8)                
                if (histoname == "effvsxy"):                   
                    pline = ROOT.TPolyLine(5,xp,yp)
                    pline.SetLineColor(2)
                    pline.SetLineWidth(2)
                    pline.Draw()
                c1.SaveAs("Run_"+runNumber+"/"+histoname+"_Run"+runNumber+".pdf")  
 
            if (histoname == "effvst3"):
                myFit = myhisto.Fit("pol0","QS")            
                print "eff media per Run ", runNumber," = ",round(myFit.Parameter(0),4), " +/- ", round(myFit.ParError(0),4)
            if (histoname == "linq"):
                #myhisto.GetXaxis().SetRangeUser(0.,20.)
                myhisto.Draw()
                c1.SaveAs("Run_"+runNumber+"/"+histoname+"_Run"+runNumber+".pdf")                
        except:
            print histoname, "for run ",runNumber," is missing"
    eff = evenEntries   / oddEntries
    ratioNumbers = round(eff,2)
    error2 = np.power(evenEntries,2)*oddEntries/np.power(oddEntries,4) + evenEntries/np.power(oddEntries,2)
    error = np.sqrt(error2)
    ratioNumbersError = round(error,2)
    #print "Run ",runNumber," Ratio between # clusters in even / odd rows ", ratioNumbers                
    return ratioNumbers, ratioNumbersError


#main part
ROOT.gROOT.SetBatch()
ROOT.gErrorIgnoreLevel = ROOT.kFatal
ROOT.gStyle.SetOptStat(0)                                                                                                                                           
ratios = []
ratiosError = []
runlist = [37692,38616,38617,38615,38614,38613,38619]
vthr = [1000,1151,1151,1570,2158,2170,3514]

#runlist = [37692]
#vthr = [1000]
for runNumber in runlist:
    mydir = "Run_"+str(runNumber)
    if os.path.exists(mydir):
        try:
            os.system('rm -rf %s' % mydir)
        except OSError:
            print ("Removal of  of the directory %s failed" % mydir)
        try:
            os.mkdir(mydir)
        except OSError:
            print ("Creation of the directory %s failed" % path)
    else:
        try:
            os.mkdir(mydir)
        except OSError:
            print ("Creation of the directory %s failed" % path)
    ratio, ratioErr = plotHisto(runNumber)
    #print runNumber, ratio, ratioErr
    ratios.append(ratio)
    ratiosError.append(ratioErr)

fig =plt.figure(figsize=(7,5), dpi=100)
ax = fig.add_subplot(111)
ax.xaxis.grid(True, which="minor")
ax.yaxis.grid(True, which="major")

plt.xlim(800,3600)
plt.ylim(0,9)
plt.xticks(fontsize = 8)
y_axis_label = "# Events cluster Row Even/Odd"
plt.ylabel(y_axis_label)
plt.xlabel("VThreshold values [e]")
plt.plot(vthr[1:],ratios[1:],'bo',markersize=6,linewidth=0)
plt.plot(vthr[0],ratios[0],'ro',markersize=6,linewidth=0)
ax.patch.set_facecolor("w")
fig.patch.set_facecolor("w")
#adding the legend
from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color="b", lw=3),
                Line2D([0], [0], color="r", lw=3)]

#ax.legend(custom_lines, ['Italy', other_country, 'France'], loc='upper left')
ax.legend(custom_lines, ['Planar', '3D'], loc='upper right')
plt.show()
