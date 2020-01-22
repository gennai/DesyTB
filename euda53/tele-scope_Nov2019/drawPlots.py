import ROOT
import array
import math
import os

mpvDict ={}

def plotHisto(runNumber):
    histolist = [
    "linqx",
    "effvsxmym",
    "effvsxmym",
    "effvsxy",
    "dutpxqvsxy",
    "linnpxvsxmym",
    "linqxvsxmym",
    "linqxvsxmymaverage",
    "linqxvsxy",
    "linq",
    "dutdxdy",
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
    c1.SetRightMargin(0.2423698)
    c1.SetFrameBorderMode(0)
    c1.SetFrameBorderMode(0) 

    meanFromMoyalHisto = 0          
    for histoname in histolist:
        myhisto = myfile.Get(histoname)
        try:
            if (histoname == "linqx"):
                    meanFromMoyalHisto = myhisto.GetMean()
                    meanFromMoyalHisto = -math.log(meanFromMoyalHisto)

            if (histoname != "effvst3" and histoname != "linq" and histoname != "linqx"):                                   
                """
                if ("linqxvsxmymaverage" in histoname):
                    px = myhisto.ProfileX()
                    px.SetTitle("Lin cluster ToT vs xmod")
                    px.GetXaxis().SetTitle("x track mod [100 #mum]")
                    px.GetXaxis().SetTitleOffset(1.2)
                    px.GetYaxis().SetTitle("ToT")                
                    px.Draw()
                    c1.SaveAs("Run_"+runNumber+"/"+histoname+"Profile_Run"+runNumber+".pdf")    
                """
                myhisto.Draw("colz")
                myhisto.GetZaxis().SetTitleOffset(1.8)
                if (histoname == "effvsxy"):                   
                    pline = ROOT.TPolyLine(5,xp,yp)
                    pline.SetLineColor(2)
                    pline.SetLineWidth(2)
                    pline.Draw()
                c1.SaveAs("Run_"+runNumber+"/"+histoname+"_Run"+runNumber+".pdf")  
            """
            if (histoname == "linnpxvsxmym"):                
                px = myhisto.ProfileX()
                px.SetTitle("Lin cluster size vs xmod")
                px.GetXaxis().SetTitle("x track mod [100 #mum]")
                px.GetXaxis().SetTitleOffset(1.2)
                px.GetYaxis().SetTitle("Cluster size")                
                px.Draw()
                c1.SaveAs("Run_"+runNumber+"/"+histoname+"Profile_Run"+runNumber+".pdf")  
            """     
            if (histoname == "effvst3"):
                myFit = myhisto.Fit("pol0","QS")            
                print "eff media per Run ", runNumber," = ",round(myFit.Parameter(0),4), " +/- ", round(myFit.ParError(0),4)
            if (histoname == "linq"):
                myhisto.GetXaxis().SetRangeUser(0.,30.)
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

                #myFitL.Print()
                meanFromMoyalHisto = meanFromMoyalHisto*myFitL.Parameter(2)
                print "MPV from fit for Run ",runNumber," = ",round(myFitL.Parameter(1),1), " and sigma = ", round(myFitL.Parameter(2),2) 
                print "Mean value for Run ",runNumber," = ",round(mean,1)
                print "MPV from Moyal ",runNumber, " = ", round(meanFromMoyalHisto,1)
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
                c1.SaveAs("Run_"+runNumber+"/"+histoname+"_Run"+runNumber+".pdf")
                mpvDict[runNumber] = [round(mean,1), round(meanFromMoyalHisto,1),round(myFitL.Parameter(1),1)]
        except:
            print histoname, "for run ",runNumber," is missing"


#main part
ROOT.gROOT.SetBatch()
ROOT.gErrorIgnoreLevel = ROOT.kFatal
ROOT.gStyle.SetOptStat(0)                                                                                                                                           

#runlist = [37673,37674,37676,37677,37691,37692,37631,37722, 37723, 37724]
runlist = [37692,37676,37674,37722,37724,37631,37683,37643]
#runlist = [37683,37643]
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

meanArray = array.array( 'd',[mpvDict[str(key)][0] for key in runlist])                                                                                                        
moyalArray = array.array( 'd',[mpvDict[str(key)][1] for key in runlist])                                                                                                        
mpvArray = array.array( 'd',[mpvDict[str(key)][2] for key in runlist])                                                                                                        
meanArrayErrors = array.array("d",[0.05 for a in range(len(runlist))])   
#print zoneArray
#print meanArray

#Plotting the various Landau approx
c2 = ROOT.TCanvas("c2","",900,600)
c2.SetFillColor(0)
c2.SetBorderMode(0)
c2.SetBorderSize(2)
c2.SetRightMargin(0.2423698)
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
gr.SetMarkerColor(4)                                                                                                                                  
gr.SetMarkerSize(markerSize)                                                                                                                                  
gr.Draw("Psame")            
gr1 = ROOT.TGraphErrors( int(len(runlist)),zoneArray,moyalArray,zoneArrayError,meanArrayErrors )                                                          
gr1.SetMarkerStyle(22)                                                                                                                                 
gr1.SetMarkerColor(2)                                                                                                                                  
gr1.SetMarkerSize(markerSize)                                                                                                                                  
gr1.Draw("Psame")            
gr2 = ROOT.TGraphErrors( int(len(runlist)),zoneArray,mpvArray,zoneArrayError,meanArrayErrors )                                                          
gr2.SetMarkerStyle(23)                                                                                                                                 
gr2.SetMarkerColor(3)                                                                                                                                  
gr2.SetMarkerSize(markerSize)                                                                                                                                  
gr2.Draw("Psame") 
xl1=.8                                                                                                                                               
yl1=0.73                                                                                                                                              
xl2=xl1+.25                                                                                                                                           
yl2=yl1+.15                                                                                                                                           
x3l1 = xl1 + 0.1                                                                                                                                      
leg1 = ROOT.TLegend(xl1,yl1,xl2,yl2)                                                                                                                  
leg1.SetBorderSize(0)   
leg1.AddEntry(gr,"Mean","p")                              
leg1.AddEntry(gr2,"MPV","p")     
leg1.AddEntry(gr1,"Moyal","p")                              
                         
leg1.Draw()

c2.SaveAs("Landau_MPVs.pdf") 