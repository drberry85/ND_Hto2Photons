print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
gROOT.ProcessLine(".L effSigma.C")
from ROOT import effSigma
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

print "Setting Initial Parameters."
HistogramNames=[["PixelBarrelSuperdZRes","TIBSuperdZRes","TOBSuperdZRes","PixelFwdSuperdZRes","TIDSuperdZRes","TECSuperdZRes"]]
HistogramNames.append(["PixelBarrelConvdZRes","TIBConvdZRes","TOBConvdZRes","PixelFwdConvdZRes","TIDConvdZRes","TECConvdZRes"])
filename="Vertex_Background.root"
pwd = "../"
file=TFile(pwd+filename)

for histlist in HistogramNames:
	names=""
	results=""
	resultsRMS=""
	yeilds=""
	for histname in histlist:
		hist=file.Get(histname)
		RMS=hist.GetRMS()
		sigma = effSigma(hist)
		results+=histname+": "+str(round(sigma,3))[:5]+" "
		resultsRMS+=histname+": "+str(round(RMS,3))[:5]+" "
		yeilds+=histname+": "+str(hist.Integral(0,hist.GetNbinsX()+1))+" "
	print "Effective Sigma: "+results
	#print "RMS: "+resultsRMS
	#print "Yeilds: "+yeilds

print "Done!"
