import ROOT
import array
import math
import os
import numpy as np

mpvDict ={}

def fitLandauGauss( x, par ):

    nn = 0
    nn=nn+1
    xbin = 1
    b1 = 0
    if( nn == 1 ):
        b1 = x[0]
   
    if( nn == 2 ):
        xbin = x[0] - b1 # bin width needed for normalization
    

    # Landau:

    invsq2pi = 0.3989422804014   #// (2 pi)^(-1/2)
    mpshift  = -0.152278298       #// Landau maximum location

    #MP shift correction:

    mpc = par[0] - mpshift * par[1] #//most probable value (peak pos)

    #//Fit parameters:
    #//par[0] = Most Probable (MP, location) parameter of Landau density
    #//par[1] = Width (scale) parameter of Landau density
    #//par[2] = Total area (integral -inf to inf, normalization constant)
    #//par[3] = Gaussian smearing

    #// Control constants
    np = 100.0      #// number of convolution steps
    sc =   5.0      #// convolution extends to +-sc Gaussian sigmas

    #// Range of convolution integral
    xlow = x[0] - sc * par[3]
    if( xlow < 0 ): xlow = 0
    xupp = x[0] + sc * par[3]

    step = (xupp-xlow) / np

    # // Convolution integral of Landau and Gaussian by sum

    sum = 0
    xx = 0
    fland = 0

    for i in range(int(np/2)):
     
        xx = xlow + ( i - 0.5 ) * step
        fland = ROOT.TMath.Landau( xx, mpc, par[1] ) / par[1]
        sum += fland * ROOT.TMath.Gaus( x[0], xx, par[3] )

        xx = xupp - ( i - 0.5 ) * step
        fland = ROOT.TMath.Landau( xx, mpc, par[1] ) / par[1]
        sum += fland * ROOT.TMath.Gaus( x[0], xx, par[3] )
  

    return par[2] * invsq2pi * xbin * step * sum / par[3] 


def f_moyal(x,par):    
    nn = 0
    dx = 1
    x0 = 0
    if( nn == 0 ):
        x0 = x[0]
    if( nn == 1 ):
        dx = x[0] - x0; # bin width needed for normalization
  
    nn=nn+1

    #Moyal: a1 * exp{-a2 ( t + a3* exp[-a4* t ] },

    t = ( x[0] - par[0] ) / par[1]
    f = np.exp( -par[2]* t - par[3]*np.exp( -par[4] * t ) )

    return par[5]*dx*0.4*f/par[1] #0.4 = 1/sqrt(2pi) norm

def plotHisto(runNumber):
    histolist = [
    "linqx",
    "linqxmoyal",
    "effvsxmym",
    "effvsxy",
    "dutpxqvsxy",
    "linnpxvsxmym",
    "linqxvsxmym",
    "linqxvsxmymaverage",
    "linqxvsxy",
    "linq",
    "dutdxdy",
    "dutdx",
    "dutdyc2",
    "effvst3",
    ]   
    if (runNumber == 37692 or runNumber == 37724): #highstat runs
        histolist.append("effvsxmymhighstat")
        histolist.append("linnpxvsxmymhighstat")
        histolist.append("linqxvsxmymhighstat")
    if (runNumber == 37692 ): #highstat runs 25x100
        histolist.append("linnpxvsxmymhighstatlargecell")

    x = [-3.5,-3.5,3.2,3.2,-3.5]
    y = [-4.7,4.7,4.7,-4.7,-4.7]

    if ( 37721 < runNumber <37725):
        x = [-3.5,-3.5,3.1,3.1,-3.5]
        y = [-4.0,4.7,4.7,-4.0,-4.0]
    
    if ( runNumber == 37631):
        x = [-3.1,-3.1,3.1,3.1,-3.1]
        y = [-4.7,4.7,4.7,-4.7,-4.7]

    if (37672 < runNumber < 37693 ):    
        x = [-4.7,-4.7,1.,1.,-4.7]
        y = [-3.5,3.1,3.1,-3.5,-3.5]


    xp = xbins = array.array('d',x)
    yp = xbins = array.array('d',y)
    
    if (runNumber != 37692):
        runNumber = str(runNumber) 
        fileToOpen = "scopeRD"+runNumber+".root"
    else:
        runNumber = str(runNumber) 
        #fileToOpen = "scopeRD"+runNumber+"_15k.root"
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
            if (histoname == "linqx"):
                    meanFromMoyalHisto = myhisto.GetMean()
                    meanFromMoyalHisto = -math.log(meanFromMoyalHisto)
                    print "mean from moyal", meanFromMoyalHisto
            if (histoname == "linqxmoyal"):
                    meanFromMoyalHistoWithMoyalWidth = myhisto.GetMean()
                    meanFromMoyalHistoWithMoyalWidth = -math.log(meanFromMoyalHistoWithMoyalWidth)
                    print "mean from moyal", meanFromMoyalHistoWithMoyalWidth

            if (histoname != "effvst3" and histoname != "linq" and histoname != "linqx" and histoname != "linqxmoyal"):                                                  
                myhisto.Draw("colz")
                myhisto.GetZaxis().SetTitleOffset(1.8)
                #if (histoname == "linqxvsxmym"):
                    #myhisto.SetMaximum(11.)
                    #myhisto.SetMinimum(9.)
                #if (histoname == "effvsxy"):                   
                #    pline = ROOT.TPolyLine(5,xp,yp)
                #    pline.SetLineColor(2)
                #    pline.SetLineWidth(2)
                #    pline.Draw()
                c1.SaveAs("Run_"+runNumber+"/"+histoname+"_Run"+runNumber+".png")  
 
            if (histoname == "effvst3"):
                myFit = myhisto.Fit("pol0","QS")            
                print "eff media per Run ", runNumber," = ",round(myFit.Parameter(0),4), " +/- ", round(myFit.ParError(0),4)
            if (histoname == "linq"):
                #myhisto.GetXaxis().SetRangeUser(0.,20.)
                myhisto.Draw()
                mean = myhisto.GetMean()
                fitxmin = 7
                fitxmax = 15
                if (runNumber == "37683"):
                    fitxmin = 9.
                    fitxmax = 15.
                if (runNumber == "37643"):
                    fitxmin = 10.
                    fitxmax = 16.
                myFitL = myhisto.Fit("landau","QS","",fitxmin,fitxmax)                
                meanFromMoyalHisto = meanFromMoyalHisto*myFitL.Parameter(2) 
                             
                ipk = myhisto.GetMaximumBin()
                aa = myhisto.GetSumOfWeights()
                xpk = myhisto.GetBinCenter(ipk)
                if runNumber == "37724" or runNumber == "37631": xpk = 11
                sm = xpk / 7


                x0 = xpk - 2.0*sm
                x9 = xpk + 3.5*sm
                
                fitF = ROOT.TF1("fitF",f_moyal,x0,x9,6)
                fitF.SetParameter( 0, xpk )
                fitF.SetParameter( 1, sm )
                fitF.SetParameter( 2, 0.5 )
                fitF.SetParameter( 3, 0.5 )
                fitF.SetParameter( 4, 0.9 )
                fitF.SetParameter( 5, aa )
                myFitMoyal = myhisto.Fit(fitF,"QSRL+","ep") 
                meanFromMoyalHistoWithMoyalWidth =  meanFromMoyalHistoWithMoyalWidth*myFitMoyal.Parameter(1)   
                fitF.SetLineColor(3)
                fitF.Draw("same")
                
                x9 = xpk + 6.0*sm
                #x9 = xpk + 3.5*sm
                if runNumber == "37724" : x9 = 17

                ns = sm
                if x0 < 0.5 : x0 = 0.5                
                fitLG = ROOT.TF1("fitLG",fitLandauGauss,x0,x9,4)
                fitLG.SetParameter( 0, xpk )
                fitLG.SetParameter( 1, sm )
                fitLG.SetParameter( 2, aa )
                fitLG.SetParameter( 3, ns ) 
                
                myFitLG = myhisto.Fit(fitLG,"QSRL+","ep")
                fitLG.SetLineColor(1)
                fitLG.Draw("same")

                print "###########"
                print "Sigma from Landau, Landau+Gaus and Moyal Fit ", runNumber, " = ", round(myFitL.Parameter(2),2), round(myFitLG.Parameter(3),2), round(myFitMoyal.Parameter(1),2)
                print "###########"
                
                print "MPV from fit for Run ",runNumber," = ",round(myFitL.Parameter(1),1), " and sigma = ", round(myFitL.Parameter(2),2) 
                print "MPV from Landau+Gaussian fit for Run ",runNumber," = ", round(myFitLG.Parameter(0),1)
                print "Mean value for Run ",runNumber," = ",round(mean,1)
                print "MPV from Moyal using Landau width ",runNumber, " = ", round(meanFromMoyalHisto,1)
                print "MPV from Moyal using Moyal width ", runNumber, " = ", round(meanFromMoyalHistoWithMoyalWidth,1)
                """
                xlabel = ROOT.TText()
                xlabel.SetNDC()
                xlabel.SetTextFont(2)
                xlabel.SetTextColor(1)
                xlabel.SetTextSize(0.04)                
                xlabel.SetTextAngle(0)
                xlabel.DrawText(0.6, 0.7, "Mean = "+str(round(mean,1)))
                xlabel.DrawText(0.6, 0.6, "Moyal= "+str(round(meanFromMoyalHisto,1)))
                xlabel.DrawText(0.6, 0.5, "MPV = "+str(round(myFitL.Parameter(1),1)))
                xlabel.DrawText(0.6, 0.4, "Sigma = "+str(round(myFitL.Parameter(2),2)))
                """
                c1.SaveAs("Run_"+runNumber+"/"+histoname+"_Run"+runNumber+".png")
                mpvDict[runNumber] = [round(mean,1), round(meanFromMoyalHisto,1),round(myFitL.Parameter(1),1), round(meanFromMoyalHistoWithMoyalWidth,1),round(myFitLG.Parameter(0),1)]
        except:
            print histoname, "for run ",runNumber," is missing"


#main part
ROOT.gROOT.SetBatch()
ROOT.gErrorIgnoreLevel = ROOT.kFatal
ROOT.gStyle.SetOptStat(0)                                                                                                                                           

#runlist = [37673,37674,37676,37677,37691,37692,37631,37722, 37723, 37724]
#runlist = [37692,37676,37674,37722,37724,37631,37683,37643]
runlist = [38808,38814,38815,38828,38829]

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

zoneArray = array.array("d",runlist)                                                                                                                  
zoneArrayError = array.array("d",[0.01 for a in range(len(runlist))])                                                                                           
deltaMoyalWMoyal = [(mpvDict[str(key)][3] - mpvDict[str(key)][4])/mpvDict[str(key)][4] for key in runlist]
deltaMoyalWLandau = [(mpvDict[str(key)][1] - mpvDict[str(key)][4])/mpvDict[str(key)][4] for key in runlist]
deltaLvsLG = [(mpvDict[str(key)][2] - mpvDict[str(key)][4])/mpvDict[str(key)][4] for key in runlist]
deltaMean = [(mpvDict[str(key)][0] - mpvDict[str(key)][4])/mpvDict[str(key)][4] for key in runlist]

deltaDistribMWM = ROOT.TH1F("deltaDistribMWM","",10,-0.15,0.15)
deltaDistribMWL = ROOT.TH1F("deltaDistribMWL","",10,-0.15,0.15)
deltaDistribMean = ROOT.TH1F("deltaDistribMean","",10,-0.15,0.15)
deltaDistribLvsLG = ROOT.TH1F("deltaDistribLvsLG","",10,-0.15,0.15)
for i in range(len(runlist)):
    deltaDistribLvsLG.Fill(deltaLvsLG[i])
    deltaDistribMWM.Fill(deltaMoyalWMoyal[i])
    deltaDistribMWL.Fill(deltaMoyalWLandau[i])
    deltaDistribMean.Fill(deltaMean[i])

c3 =  ROOT.TCanvas("c2","",900,600)
deltaDistribLvsLG.GetXaxis().SetTitle("(Estimation - MPV)/MPV")
deltaDistribLvsLG.SetLineColor(3)
deltaDistribLvsLG.SetFillColor(3)
deltaDistribLvsLG.Draw()
deltaDistribMean.SetLineColor(1)
deltaDistribMean.SetFillColor(1)
deltaDistribMean.SetFillStyle(22)
deltaDistribMean.Draw("same")
deltaDistribMWM.SetLineColor(6)
deltaDistribMWM.SetFillColor(6)
deltaDistribMWM.Draw("same")
deltaDistribMWL.SetLineColor(2)
deltaDistribMWL.SetFillColor(2)
deltaDistribMWL.Draw("same")
xl1=.15                                                                                                                                               
yl1=0.73                                                                                                                                              
xl2=xl1+.25                                                                                                                                           
yl2=yl1+.15                                                                                                                                           
x3l1 = xl1 + 0.1                                                                                                                                      
leg2 = ROOT.TLegend(xl1,yl1,xl2,yl2)                                                                                                                  
leg2.SetBorderSize(0)   
leg2.AddEntry(deltaDistribMean,"Mean","l")       
leg2.AddEntry(deltaDistribLvsLG,"Simple Landau","l")                         
leg2.AddEntry(deltaDistribMWL,"Moyal w/ Landau #sigma","l") 
leg2.AddEntry(deltaDistribMWM,"Moyal w/ Moyal #sigma","l")   
                            
                         
leg2.Draw()

c3.SaveAs("DeltaDistrib.png")   

meanArray = array.array( 'd',[mpvDict[str(key)][0] for key in runlist])                                                                                                        
moyalArray = array.array( 'd',[mpvDict[str(key)][1] for key in runlist])                                                                                                        
mpvArray = array.array( 'd',[mpvDict[str(key)][2] for key in runlist])                                                                                                        
mpvLGArray = array.array( 'd',[mpvDict[str(key)][4] for key in runlist])                                                                                                        
moyalWithMoyalSigmaArray = array.array( 'd',[mpvDict[str(key)][3] for key in runlist])                                                                                                        
meanArrayErrors = array.array("d",[0.05 for a in range(len(runlist))])   


#Plotting the various Landau approx
c2 = ROOT.TCanvas("c2","",900,600)
c2.SetFillColor(0)
c2.SetBorderMode(0)
c2.SetBorderSize(2)
c2.SetRightMargin(0.15423698)
c2.SetFrameBorderMode(0)
c2.SetFrameBorderMode(0) 
h1 = ROOT.TH2F("h1","",len(runlist),37620,37740,100,8.0,15.0)                                                                                                          
h1.GetXaxis().SetTitle("Runs")               
h1.GetXaxis().SetTitleOffset(1.2)                                                                            
h1.GetXaxis().SetLabelOffset(0.01)                                                                            
h1.GetYaxis().SetTitleOffset(1.2)

markerSize = 1.3
h1.GetYaxis().SetTitle("Mean/MPV/Moyal")                                                                                                        
h1.Draw()                                                                                                                                             
gr = ROOT.TGraphErrors( int(len(runlist)),zoneArray,meanArray,zoneArrayError,meanArrayErrors )                                                          
gr.SetMarkerStyle(21)                                                                                                                                 
gr.SetMarkerColor(1)  
gr.SetLineColor(1)                                                                                                                                       
gr.SetMarkerSize(markerSize)                                                                                                                                  
gr.Draw("lPSame")            
gr1 = ROOT.TGraphErrors( int(len(runlist)),zoneArray,moyalArray,zoneArrayError,meanArrayErrors )                                                          
gr1.SetMarkerStyle(22)                                                                                                                                 
gr1.SetMarkerColor(2)   
gr1.SetLineColor(2)                                                                                                                                      
gr1.SetMarkerSize(markerSize)                                                                                                                                  
gr1.Draw("lPSame")            
gr2 = ROOT.TGraphErrors( int(len(runlist)),zoneArray,mpvArray,zoneArrayError,meanArrayErrors )                                                          
gr2.SetMarkerStyle(23)                                                                                                                                 
gr2.SetMarkerColor(3)  
gr2.SetLineColor(3)                                                                                                                                
gr2.SetMarkerSize(markerSize)                                                                                                                                  
gr2.Draw("lPSame") 
gr3 = ROOT.TGraphErrors( int(len(runlist)),zoneArray,moyalWithMoyalSigmaArray,zoneArrayError,meanArrayErrors )                                                          
gr3.SetMarkerStyle(24)                                                                                                                                 
gr3.SetMarkerColor(6)  
gr3.SetLineColor(6)                                                                                                                                
gr3.SetMarkerSize(markerSize)                                                                                                                                  
gr3.Draw("lPSame") 
gr5 = ROOT.TGraphErrors( int(len(runlist)),zoneArray,mpvLGArray,zoneArrayError,meanArrayErrors )                                                          
gr5.SetMarkerStyle(26)                                                                                                                                 
gr5.SetMarkerColor(7)  
gr5.SetLineColor(7)                                                                                                                                
gr5.SetMarkerSize(markerSize)                                                                                                                                  
gr5.Draw("lPSame") 
gr3 = ROOT.TGraphErrors( int(len(runlist)),zoneArray,moyalWithMoyalSigmaArray,zoneArrayError,meanArrayErrors )                                                          
gr3.SetMarkerStyle(24)                                                                                                                                 
gr3.SetMarkerColor(6)  
gr3.SetLineColor(6)                                                                                                                                
gr3.SetMarkerSize(markerSize)                                                                                                                                  
gr3.Draw("lPSame") 
xl1=.77                                                                                                                                               
yl1=0.73                                                                                                                                              
xl2=xl1+.25                                                                                                                                           
yl2=yl1+.15                                                                                                                                           
x3l1 = xl1 + 0.1                                                                                                                                      
leg1 = ROOT.TLegend(xl1,yl1,xl2,yl2)                                                                                                                  
leg1.SetBorderSize(0)   
leg1.AddEntry(gr,"Mean","p")                              
leg1.AddEntry(gr2,"MPV Landau","p")     
leg1.AddEntry(gr5,"MPV L+G","p")     
leg1.AddEntry(gr1,"Moyal w/ Landau #sigma","p") 
leg1.AddEntry(gr3,"Moyal w/ Moyal #sigma","p")                              
                         
leg1.Draw()

c2.SaveAs("Landau_MPVs.png") 