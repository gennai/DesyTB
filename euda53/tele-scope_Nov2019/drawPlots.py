import ROOT
import array



def plotHisto(runNumber):
    histolist = [
    "effvsxmym",
    "effvsxmym",
    "effvsxy",
    "dutpxqvsxy",
    "linnpxvsxmym",
    "linqxvsxmym",
    "linqxvsxmymaverage",
    "linqxvsxy",
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
          
    for histoname in histolist:
        myhisto = myfile.Get(histoname)
        try:
            if (histoname != "effvst3"):    
                #if ("linqxvsxmym" in histoname ):
                #    myhisto.SetMinimum(4.5)
                #    myhisto.SetMaximum(12.0)
                #if ("linqxvsxmymaverage" in histoname):
                #    myhisto.SetMinimum(8.5)
                myhisto.Draw("colz")
                myhisto.GetZaxis().SetTitleOffset(1.8)
                if (histoname == "effvsxy"):                   
                    pline = ROOT.TPolyLine(5,xp,yp)
                    pline.SetLineColor(2)
                    pline.SetLineWidth(2)
                    pline.Draw()
                c1.SaveAs(histoname+"_Run"+runNumber+".pdf")  
            if (histoname == "linnpxvsxmym"):                
                px = myhisto.ProfileX()
                px.SetTitle("Lin cluster size vs xmod")
                px.GetXaxis().SetTitle("x track mod [100 #mum]")
                px.GetXaxis().SetTitleOffset(1.2)
                px.GetYaxis().SetTitle("Cluster size")                
                px.Draw()
                c1.SaveAs(histoname+"Profile_Run"+runNumber+".pdf")    
            if (histoname == "effvst3"):
                myFit = myhisto.Fit("pol0","QS")
                print "eff media per Run ", runNumber," = ",round(myFit.Parameter(0),4), " +/- ", round(myFit.ParError(0),4)
        except:
            print histoname, "for run ",runNumber," is missing"


#main part
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)                                                                                                                                           

#runlist = [37673,37674,37676,37677,37691,37692,37631,37722, 37723, 37724]
runlist = [37692,37676,37674,37722,37724,37631]
#runlist = [37565]
for runNumber in runlist:
    plotHisto(runNumber)
