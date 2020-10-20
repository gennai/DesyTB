import ROOT
import numpy as np
import matplotlib.pyplot as plt
import sys
from math import sqrt

ROOT.gROOT.SetBatch()
histoname_ref = "dutdyvsxm"
extension = ""
try:
    sys.argv[1]
except:
    print "histoname = ", histoname_ref
else:
    histoname_ref += str(sys.argv[1])
    extension = str(sys.argv[1])
    print "histoname = ", histoname_ref

vsbias = True

runlist = [39313,39330,39324,39335,39334,39333,39321,39322,39323]
vThr = [1150,1150,1050,1400,1150,950,1315,1500,1228]
vbias = [200,800,500,500,500,500,500,500,500]
if vsbias:
    runlist= [38868,39311,39313,39316,39318,39320,39326,39328,39330]
    vbias = [70,100,200,300,400,500,600,700,800]
    vThr =  vbias
ratios = []
ratios_err = []
ratios_fit = []
ratios_fit_err = []
deltas = []
deltas_err = []
c1 = ROOT.TCanvas("c1", "", 800,600)
for (thr, run, bias) in zip(vThr,runlist, vbias):
    fileToOpen = "scopeRD"
    fileToOpen += str(run)+".root"
    myfile = ROOT.TFile(fileToOpen)
    if run == 38868:
        histoname = "dutdxvsym"
    else:
        histoname = histoname_ref

    myhistobond =   myfile.Get(histoname+"bonding")
    myhistonobond = myfile.Get(histoname+"nobonding")
    myhistobond.Draw()
    f2 = ROOT.TF1("f2","[1]*sin(x*6.28/100+[0])+[2]")
    myhistobond.Fit("f2","QS","",0,100)
    #c1.SaveAs("fitBonding_"+extension+"_run"+str(run)+".pdf")
    delta_fit = abs(f2.GetParameter(1))*2
    nbins = myhistobond.GetNbinsX()
    hrms = ROOT.TH1F("hrms","",100,-0.01, 0.01)
    for i in range(nbins):
        hrms.Fill(myhistonobond.GetBinContent(i))
    hrms.Draw()
    print run, myhistobond.GetEntries(), hrms.GetEntries()
    f1 = hrms.Fit("gaus","QS","")
    #print "rms from fit ", f1.Parameter(2)
    #c1.SaveAs("fitNoBonding_"+extension+"_run"+str(run)+".pdf")
    err1 = myhistobond.GetBinError(myhistobond.GetMaximumBin())
    err2 = myhistobond.GetBinError(myhistobond.GetMinimumBin())
    err3 = myhistonobond.GetBinError(myhistonobond.GetMaximumBin())
    err4 = myhistonobond.GetBinError(myhistonobond.GetMinimumBin())
    num = myhistobond.GetMaximum() - myhistobond.GetMinimum()
    den = myhistonobond.GetMaximum() - myhistonobond.GetMinimum()
    ratio = num/den
    ratios.append(ratio)
    fit_rms = round(f1.Parameter(2),4)
    error = err1*sqrt(2)/den
    ratios_err.append(error)
    deltas.append(delta_fit*1E3)
    deltas_err.append(f2.GetParError(1)*1E3)
    ratios_fit.append(delta_fit/(2*fit_rms))
    error_fit = f2.GetParError(1)*f2.GetParError(1)/(delta_fit*delta_fit)
    error_fit += f1.ParError(2)*f1.ParError(2)/( fit_rms*fit_rms )
    error_fit = sqrt(error_fit)
    error_fit = error_fit*delta_fit/(2*fit_rms)
    ratios_fit_err.append(error_fit)
    print run, thr, round(delta_fit,6), fit_rms, round(delta_fit/(2*fit_rms),4), "+/-", round(error_fit,4)
    """
    print " Run ",run," Vthr ", thr, " Vbias ",bias," ratio ", ratio, \
          " delta bond ", round(myhistobond.GetMaximum() - myhistobond.GetMinimum(),4), \
          " delta fit ", round(delta_fit, 6) ,\
          " delta no bond ", round(myhistonobond.GetMaximum() - myhistonobond.GetMinimum(),4), \
          " rms no bond ", round(sigma_nobond,6)
    """
    #print "Run ", run, " ratio ", ratio, " fit ", delta_fit/hrms.GetRMS()
fig =plt.figure(figsize=(7,5), dpi=100)
ax = fig.add_subplot(111)
ax.xaxis.grid(True, which="minor")
ax.yaxis.grid(True, which="major")
plt.xlim(900,1600)
plt.ylim(0.,6.)
plt.xticks(fontsize = 8)
y_axis_label = " Ratio = Delta Bond / Delta NoBond"
plt.ylabel(y_axis_label)
plt.xlabel("Threshold [electrons]")
#plt.errorbar(vThr,ratios,yerr=ratios_err,fmt="bo",markersize=6,linewidth=0,elinewidth=1 )
#plt.errorbar(vThr[0],ratios[0],yerr=ratios_err[0],fmt="go",markersize=6,linewidth=0,elinewidth=2 )
#plt.errorbar(vThr[1],ratios[1],yerr=ratios_err[1],fmt="ro",markersize=6,linewidth=0,elinewidth=2 )
ax.patch.set_facecolor("w")
fig.patch.set_facecolor("w")
#adding the legend
#plt.show()
from matplotlib.lines import Line2D

if vsbias == True:
    plt.xlim(50,850)
    plt.subplot(ax)
    plt.xlabel("V bias [V]")
    plt.errorbar(vbias[0:],ratios_fit[0:],yerr=ratios_fit_err[0:],fmt='bs',markersize=6,linewidth=0,elinewidth=1, label="irrad" )
    plt.errorbar(vbias[0],ratios_fit[0],yerr=ratios_fit_err[0],fmt='rs',markersize=6,linewidth=0,elinewidth=1, label = "new" )
    plt.legend(loc='best',numpoints=1)
#    ax1 = fig.add_subplot(212, sharex=ax, sharey=ax)
#    plt.subplot(ax1)
#    plt.errorbar(vbias[0:],deltas[0:],yerr=deltas_err[0:],fmt='bs',markersize=6,linewidth=0,elinewidth=1, label="irrad" )
#    plt.errorbar(vbias[0],deltas[0],yerr=deltas_err[0],fmt='rs',markersize=6,linewidth=0,elinewidth=1, label = "new" )
    plt.savefig("ratioVsBias_"+histoname+".pdf")
else:   
    plt.errorbar(vThr,ratios_fit,yerr=ratios_fit_err,fmt='bv',markersize=6,linewidth=0,elinewidth=1 )
    plt.errorbar(vThr[0],ratios_fit[0],yerr=ratios_fit_err[0],fmt='gv',markersize=6,linewidth=0,elinewidth=1 )
    plt.errorbar(vThr[1],ratios_fit[1],yerr=ratios_fit_err[1],fmt='rv',markersize=6,linewidth=0,elinewidth=1 )
    custom_lines = [Line2D([0], [0], color="b", lw=3),
                Line2D([0], [0], color="r", lw=3),
                Line2D([0], [0], color="g", lw=3)]
    ax.legend(custom_lines, ['Vbias = 500 V', 'Vbias = 800 V','Vbias = 200 V'], loc='lower left')
    plt.savefig("ratioVsThreshold_"+histoname+".pdf")