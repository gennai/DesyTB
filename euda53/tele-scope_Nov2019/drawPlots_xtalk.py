import ROOT
import array
import math
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse



def plotHisto(runNumber):
    histolist = [
    "linqCell1",
    "linqCell2",
    "linqx",
    "linqxvsym",
    "linqxmoyal",
    "effvsxmym",
    "effvsxy",
    "dutpxqvsxy",
    "linnpxvsxmym",
    "linqxvsxmym",
    "linnpxvsxm",
    "linnpxvsxmnobonding",
    "linqxvsxmymaverage",
    "dutdxdy",
    "dutdyvsxmym",
    "dutdyvsxmym1Cell",
    "dutdyvsxm",
    "dutdyvsym",
    "dutdx",
    "dutdxvsxm",
    "dutdyc2",
    "sixdxc",
    "linq",
    "effvsxmym",
    "effvsxy",
    "linqxvsxy",
    "linqxvsx",    
    "linq",
    "linrowmin2",
    "dutdxdy",
    "linqxvsy",
#    "dutdyc2",
#    "linqMinOverqCluster",
#    "linqMinOverqClusterAlongRows",
#    "linqMinOverqClusterAlongColumns",
    "effvst3",
    "linnrow1odd",
#    "sixdtx",
#    "sixdty",
#    "sixdtyLargeClusters", #slope on x for driplets
#    "sixdtxLargeClusters", #slope on y for driplets
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
    myHistoRes = ROOT.TH1F("Resolution","DUT - track dx;DUT cluster - track #Deltax [mm];DUT clusters",500, -0.2, 0.2 )
    if runNumber == "39321" or runNumber == "38844" or runNumber == "38880" or runNumber == "38881" or runNumber == "38903" or runNumber == "38873":        
        dutdxPaired = myfile.Get("dutdxPaired")
        dutdxPairedCls1 = myfile.Get("dutdxPairedCls1")
        dutdxPairedCls2 = myfile.Get("dutdxPairedCls2")
        dutdxPairedCls3 = myfile.Get("dutdxPairedCls3")
        dutdxUnpaired = myfile.Get("dutdxUnpaired")
        dutdxUnpairedCls1 = myfile.Get("dutdxUnpairedCls1")
        dutdxUnpairedCls2 = myfile.Get("dutdxUnpairedCls2")
        dutdxUnpairedCls3 = myfile.Get("dutdxUnpairedCls3")
        
        myHistoRes.GetXaxis().SetRangeUser(-0.05,0.05)
        myHistoRes.SetMaximum(0.07)
        myHistoRes.Draw()
        dutdxUnpaired.DrawNormalized("same")
        dutdxPaired.SetLineColor(2)
        dutdxPaired.DrawNormalized("same")
        xl1=.15                                                                                                                                               
        yl1=0.73                                                                                                                                              
        xl2=xl1+.25                                                                                                                                           
        yl2=yl1+.15                                                                                                                                           
        x3l1 = xl1 + 0.1                                                                                                                                      
        leg2 = ROOT.TLegend(xl1,yl1,xl2,yl2)                                                                                                                  
        leg2.SetBorderSize(0)   
        leg2.AddEntry(dutdxUnpaired,"Unpaired rows","l")  
        leg2.AddEntry(dutdxPaired,"Paired rows","l")  
        leg2.Draw()        
        c1.SaveAs("Run_"+runNumber+"/dutdxAll_Run"+runNumber+".png")
        myHistoRes.SetMaximum(0.12)
        myHistoRes.Draw()
        dutdxUnpairedCls1.DrawNormalized("same")
        dutdxPairedCls1.SetLineColor(2)
        dutdxPairedCls1.DrawNormalized("same")
        leg2.Draw()        
        c1.SaveAs("Run_"+runNumber+"/dutdxCls1_Run"+runNumber+".png")
        myHistoRes.SetMaximum(0.07)
        myHistoRes.Draw()
        dutdxUnpairedCls2.DrawNormalized("same")
        dutdxPairedCls2.SetLineColor(2)
        dutdxPairedCls2.DrawNormalized("same")
        leg2.Draw()        
        c1.SaveAs("Run_"+runNumber+"/dutdxCls2_Run"+runNumber+".png")
        myHistoRes.Draw()
        dutdxUnpairedCls3.DrawNormalized("same")
        dutdxPairedCls3.SetLineColor(2)
        dutdxPairedCls3.DrawNormalized("same")
        leg2.Draw()        
        c1.SaveAs("Run_"+runNumber+"/dutdxCls3_Run"+runNumber+".png")
    
    for histoname in histolist:
        myhisto = myfile.Get(histoname)
        try:
           # if ("sixdt" in histoname):
            #    print histoname,"\n"
             #   print "mean and rms ", round(myhisto.GetMean(),6), round(myhisto.GetRMS(),6),"\n"
           
            if(histoname =="linnrow1odd"):
                oddEntries = myhisto.GetBinContent(2)
            if(histoname =="linnrow1eve"):
                evenEntries = myhisto.GetBinContent(2)    
            if (histoname =="linqMinOverqCluster"):
                #myhisto.Rebin(2)
                #myhisto.GetXaxis().SetRangeUser(0.04,0.1)
                myhisto.Draw()
                c1.SaveAs("Run_"+runNumber+"/"+histoname+"_Run"+runNumber+".png")   
                #print myhisto.GetMean()
                #myhisto.SetMaximum(35000)
            if (histoname =="linrowmin2"):
                myhisto.GetXaxis().SetRangeUser(100.,140.)
                myhisto.Draw()
                c1.SaveAs("Run_"+runNumber+"/"+histoname+"_Run"+runNumber+".png")              
            if (histoname != "effvst3" and histoname != "linq" and histoname != "linrowmin2" and histoname != "linqMinOverqCluster"):                                                  
                myhisto.Draw("colz")
                myhisto.GetZaxis().SetTitleOffset(1.8)
                if (histoname == "linqxvsym" and runNumber != "39313" and runNumber != "38868"):
                    myhisto.SetMinimum(10)                                        
                if (histoname == "linqxvsxmym" and runNumber != "38868"):
                    myhisto.SetMaximum(12)
                if (histoname == "effvsxy"):                   
                    pline = ROOT.TPolyLine(5,xp,yp)
                    pline.SetLineColor(2)
                    pline.SetLineWidth(2)
                    pline.Draw()
                c1.SaveAs("Run_"+runNumber+"/"+histoname+"_Run"+runNumber+".png")  
            if (histoname == "linqxvsy"):
                    myFit = myhisto.Fit("pol1","QS")  
                    print "Slope  ", runNumber," = ",round(myFit.Parameter(1),3), " +/- ", round(myFit.ParError(1),3)
                    myhisto.Draw()
                    c1.SaveAs("Run_"+runNumber+"/"+histoname+"_Run"+runNumber+".png")                               
            if (histoname == "linqxvsx"):
                    myhisto.Draw()
                    c1.SaveAs("Run_"+runNumber+"/"+histoname+"_Run"+runNumber+".png")  
            if (histoname == "effvst3"):
                myFit = myhisto.Fit("pol0","QS")            
                print "eff media per Run ", runNumber," = ",round(myFit.Parameter(0),4), " +/- ", round(myFit.ParError(0),4)                        
            if (histoname == "linq" or histoname == "linqCell1" or histoname == "linqCell2"):
                #myhisto.GetXaxis().SetRangeUser(0.,20.)
                myhisto.Draw()
                #print histoname, myhisto.GetEntries()
                c1.SaveAs("Run_"+runNumber+"/"+histoname+"_Run"+runNumber+".png")                
        except:
            print histoname, "for run ",runNumber," is missing"
    eff = evenEntries   / (oddEntries+evenEntries)
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


#runlist = [38871,38872,38873,38874,38875,38876,38877,38878,38879,38880]
#runlist = [38903,38844,38880,38881,38873]
#runlist = [38868, 39311,39313,39324,39330, 39334, 39335]
#runlist.extend(runlist2)
#vthr = [1000]*len(runlist)
#runlist = [40418,40419,40420,40421,40423,40425,40457]
#vbias = [250,350,400,450,500,600,700]
runlist = [35472,35475,35476,35478,35479,35481,35484]
vbias = [300]*len(runlist)
vthr = vbias
#vthr = [10,20,30,40,50,60,70,80,90,100]

#runlist = [37692]
#vthr = [1000]
for runNumber in runlist:
    mydir = "Run_"+str(runNumber)
    if not os.path.exists(mydir):
        try:
            os.mkdir(mydir)
        except OSError:
            print ("Creation of the directory %s failed" % path)
    ratio, ratioErr = plotHisto(runNumber)
    print ratio 
    #print runNumber, ratio, ratioErr
    ratios.append(ratio)
    ratiosError.append(ratioErr)

fig =plt.figure(figsize=(7,5), dpi=100)
ax = fig.add_subplot(111)
ax.xaxis.grid(True, which="minor")
ax.yaxis.grid(True, which="major")

plt.xlim(900,1500)
plt.ylim(0.4,1)
plt.xticks(fontsize = 8)
y_axis_label = " #Even / (#Even + #Odd) "
plt.ylabel(y_axis_label)
#plt.xlabel("Threshold [electrons]")
plt.xlabel("Bias [V]")
plt.plot(vthr[0:],ratios[0:],'bo',markersize=6,linewidth=0)
#plt.plot(vthr[1],ratios[1],'ro',markersize=6,linewidth=0)
plt.plot(vthr[0],ratios[0],'go',markersize=6,linewidth=0)
ax.patch.set_facecolor("w")
fig.patch.set_facecolor("w")
#adding the legend
from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color="b", lw=3),
                Line2D([0], [0], color="r", lw=3)]

#ax.legend(custom_lines, ['Italy', other_country, 'France'], loc='upper left')
#ax.legend(custom_lines, ['Planar', '3D'], loc='upper right')
#plt.show()
plt.savefig("xtalk_vs_vbias.pdf")
