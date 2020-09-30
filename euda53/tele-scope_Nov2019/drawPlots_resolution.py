import ROOT
import array
import math
import os
import numpy as np
import matplotlib.pyplot as plt
G = 1.

def tp0Fit( x, par ):

    nn = 0
    nn = nn+1
    dx = 0.1
    b1 = 0
    if nn == 1 : b1 = x[0]
    if nn == 2 : dx = x[0] - b1

    # Mean and width:

    xm = par[0]
    t = ( x[0] - xm ) / par[1]
    tt = t*t

    # exponent:

    rn = par[2]
    xn = 0.5 * ( rn + 1.0 )

    # Normalization needs Gamma function:

    pk = 0.0

    if rn > 0.0 and math.fabs( xn * ROOT.TMath.Log( 1.0 + tt/rn ) ) < 333:

        pi = 3.14159265358979323846
        aa = ((dx / par[1] )/ math.sqrt(rn*pi)) * ROOT.TMath.Gamma(xn) / ROOT.TMath.Gamma(0.5*rn)

        pk = G * par[3] * aa * ROOT.TMath.Exp( -xn * ROOT.TMath.Log( 1.0 + tt/rn ) )

        # lim n->inf (1+a/n)^n = e^a

    

    return pk + par[4]

irrad = True
tele = "six"
if irrad:
    tele = "trid"
def plotHisto(runNumber):
    if irrad:
        histolist = [
      
        #for irradiated
        "dutdyc",
        "tridyc"    
        ]
    else:
        histolist = [
        #For june 2020 TB
        "dutdxc",
        "sixdxc"
        ]    

    sigmaDUT = -0.5
    sigmaSIX = -0.5
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
            
            dx = myhisto.GetBinWidth(1)
            nmax = myhisto.GetBinContent(myhisto.GetMaximumBin())
            xmax = myhisto.GetBinCenter(myhisto.GetMaximumBin())
            

            nb = myhisto.GetNbinsX()
            x1 = myhisto.GetBinCenter(1)            
            x9 = myhisto.GetBinCenter(nb)
            i1 = myhisto.FindBin(x1)
            i9 = myhisto.FindBin(x9)
            n1 = myhisto.GetBinContent(i1)
            n9 = myhisto.GetBinContent(i9)
            bg = 0.5*(n1+n9)
          
               
           
            xhwhm = x1
            for ii in range(i1,i9+1):              
                if myhisto.GetBinContent(ii)-bg > 0.5*(nmax-bg) :
                   xhwhm = myhisto.GetBinCenter(ii) 
            
            
            hwhm = xhwhm-xmax          
            if hwhm < dx : 
                hwhm = dx 
         
            nn = (3.*nmax*hwhm)/dx
            nn =nn / G
         
            #print "fit from ", x1," to ", x9, " bins ", i1, " to ", i9, " = ", i9-i1+1

            #myhisto.GetXaxis().SetRangeUser(x1, x9)
            #myhisto.Draw()
            
            fitF = ROOT.TF1("fitF",tp0Fit,x1,x9,5)
            fitF.SetParName( 0, "mean" )
            fitF.SetParName( 1, "sigma" )
            fitF.SetParName( 2, "nu" )
            fitF.SetParName( 3, "area" )
            fitF.SetParName( 4, "BG" )
      
            fitF.SetParameter( 0, xmax )
            fitF.SetParameter( 1, hwhm)
            fitF.SetParameter( 2, 7.5 )
            fitF.SetParameter( 3, nn )
            fitF.SetParameter( 4, bg)

            #print "start:", " max " , nmax, " at " , xmax, ", hwhm " , hwhm, ", bg " , bg, ", area " , nn
            if "dut" in str(histoname): 
                #print histoname
                fitDUTF = ROOT.TF1("fitDUTF",tp0Fit,x1,x9,5)
                fitDUTF.SetParName( 0, "mean" )
                fitDUTF.SetParName( 1, "sigma" )
                fitDUTF.SetParName( 2, "nu" )
                fitDUTF.SetParName( 3, "area" )
                fitDUTF.SetParName( 4, "BG" )
        
                fitDUTF.SetParameter( 0, xmax )
                fitDUTF.SetParameter( 1, hwhm)
                fitDUTF.SetParameter( 2, 7.5 )
                fitDUTF.SetParameter( 3, nn )
                fitDUTF.SetParameter( 4, bg)  
                #myhisto.GetXaxis().SetRangeUser(-0.02,0.02)
                #myhisto.Draw()             
                #myFitDUT = myhisto.Fit("gaus","QSR","",-0.01,0.01)
                #sigmaDUT =myFitDUT.Parameter(2) 
                #myhisto.Sumw2()
                fitDUTF.SetParLimits(0,-0.2,0.2)
                fitDUTF.SetParLimits(2,3.5,10.)
                #fitDUTF.FixParameter(0,xmax)
                myFitDUT = myhisto.Fit("fitDUTF","QSRB","",-0.2,0.2)                
                sigmaDUT =myFitDUT.Parameter(1) 
                #print "mean ",myFitDUT.Parameter(0), " sigma ", myFitDUT.Parameter(1), " nu ", myFitDUT.Parameter(2), " area ", myFitDUT.Parameter(3), "  bg ", myFitDUT.Parameter(4)                
            if tele in str(histoname):
                fitF.SetParLimits(2,3.5,10.)
                myFit = myhisto.Fit("fitF","QSRB")                            
                sigmaSIX =myFit.Parameter(1) 
                #print "mean ",myFit.Parameter(0), " sigma ", myFit.Parameter(1), " nu ", myFit.Parameter(2), " area ", myFit.Parameter(3), "  bg ", myFit.Parameter(4)
        except:
            print histoname, "for run ",runNumber," is missing"
        c1.SaveAs("Run_"+runNumber+"/"+histoname+"_Run"+runNumber+".png")
    if sigmaDUT >0. and sigmaSIX > 0.: 
        deltaRes2 = sigmaDUT*sigmaDUT - sigmaSIX*sigmaSIX/4
        if irrad:
            deltaRes2 = sigmaDUT*sigmaDUT - sigmaSIX*sigmaSIX
        if deltaRes2 > 0:
            deltaRes = math.sqrt(deltaRes2)
            print "Run number ",runNumber," sigma DUT ", round(sigmaDUT*1000,2), " sigma telescope ",round(sigmaSIX*1000,2), " sigma true ", round(deltaRes*1000,2)
            return deltaRes
        else:
            print "ERROR on six/dut resolution!"
            return -0.0005
    else:
        return -0.0005
#main part
ROOT.gROOT.SetBatch()
ROOT.gErrorIgnoreLevel = ROOT.kFatal
ROOT.gStyle.SetOptStat(0)                                                                                                                                           



#CHIP 3D
runlist =[39321]
angle = [0]
pngfilename = "resolution_planar_irrad.png"

#w6067_BFFP_01, 25x100, FBK, Planar, Fresh
#runlist = [38881,38882,38883,38884,38885,38886,38887,38889,38890,38891,38892]
#angle = [0,5,7,8,9,10,11,12,14,16,18]
#pngfilename = "resolution_Planar_FP.png"

#w6067_BITE_02, 25x100, FBK, Planar, Fresh	
#runlist = [38844,38845,38852,38853,38854,38855,38856,38862,38863,38865,38866]
#angle = [0,4.5,6.5,8.5,9.5,10.5,11.5,12,14,16,18]
#pngfilename = "resolution_Planar_NO_FP.png"

#CHIP 3D
#runlist =[38903,38904,38905,38906,38907,38908,38909,38910,38911,38912,38914]
#angle = [0,4,6,8,9,10,11,12,14,16,18]
#pngfilename = "resolution_3D.png"

plotType = ""

#Planar VBias
#runlist = [38871,38872,38873,38874,38875,38876,38877,38878,38879,38880]
#angle = [10,20,30,40,50,60,70,80,90,100]
#pngfilename = "resolution_vs_VBias.png"
plotType = "vbias"

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

plt.xlim(min(angle)-1,max(angle)+1)    
#plt.ylim(min(sigmax)-0.5,max(sigmax)+0.5)
plt.ylim(1.0,6.0)
plt.xticks(fontsize = 8)
y_axis_label = "Resolution [um] "
plt.ylabel(y_axis_label)
plt.xlabel("Tilt angle [deg]")
if plotType == "vbias":
    plt.xlabel("V bias [V]")
plt.xticks(np.arange(0, max(angle)+1, step=10))
plt.plot(angle,sigmax,'bo',markersize=6,linewidth=0)
ax.patch.set_facecolor("w")
fig.patch.set_facecolor("w")
#adding the legend
from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color="b", lw=3),
                Line2D([0], [0], color="r", lw=3)]

#ax.legend(custom_lines, ['Italy', other_country, 'France'], loc='upper left')
#ax.legend(custom_lines, ['Planar', '3D'], loc='upper right')

plt.savefig(pngfilename)
plt.show()
