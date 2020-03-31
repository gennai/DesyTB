import ROOT
import array
import math
import os
import numpy as np


def plotHisto(runNumber):
    histolist = [
    "linq",
    "effvsxmym",
    "effvsxy",
    "dutpxqvsxy",
    "linnpxvsxmym",
    "linqxvsxmym",
    "linqxvsxmymaverage",
    "linqxvsxy",
    "linq",
    "linrowmin2",
    "dutdxdy",
    "dutdyc2",
    "linqMinOverqCluster",
    "effvst3",
    "linnrow1odd",
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

    meanFromMoyalHisto = 0 
    meanFromMoyalHistoWithMoyalWidth = 0         
    for histoname in histolist:
        myhisto = myfile.Get(histoname)
        try:
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


#main part
ROOT.gROOT.SetBatch()
ROOT.gErrorIgnoreLevel = ROOT.kFatal
ROOT.gStyle.SetOptStat(0)                                                                                                                                           

#runlist = [37692,38613,38614,38616,38617,38619]
runlis = [37692]

runlist.sort()

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

    plotHisto(runNumber)

