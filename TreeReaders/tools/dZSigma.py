"Loading Root..."

def Sigma(hist,width):
	meanbin=hist.FindBin(hist.GetMean())
	HistInt=0
	counter=0
	while HistInt<width:
		counter += 1
		lowerbin = meanbin-counter
		upperbin = meanbin+counter
		if lowerbin<0: lowerbin=0
		if upperbin>hist.GetNbinsX(): upperbin=hist.GetNbinsX()
		HistInt = hist.Integral(lowerbin,upperbin)/hist.Integral()
		#print "%i %i %i %f" %(counter,lowerbin,upperbin,HistInt)
	if len(hist.GetName())<12:
		print "%s:\t\t%f\t%0.3f\t\t%i\t\t%i" %(hist.GetName(),hist.GetBinLowEdge(upperbin+1)-hist.GetBinLowEdge(lowerbin),HistInt,lowerbin,upperbin)
	else:
		print "%s:\t%f\t%0.3f\t\t%i\t\t%i" %(hist.GetName(),hist.GetBinLowEdge(upperbin+1)-hist.GetBinLowEdge(lowerbin),HistInt,lowerbin,upperbin)

#import pdb; pdb.set_trace()
from sys import argv
from ROOT import *
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

filename = argv[1]
histlist = ["PixelBarrelNewPVdZ","TIBNewPVdZ","TOBNewPVdZ","PixelFwdNewPVdZ","TIDNewPVdZ","TECNewPVdZ","PixelBarrelLineardZ","TIBLineardZ","TOBLineardZ","PixelFwdLineardZ","TIDLineardZ","TECLineardZ"]
sigmas = [0.68,0.90]
DataFile = TFile(filename)

for sigma in sigmas:
	print "\nSigma %0.3f\t\tWidth (cm)\tActual Sigma\tUpper Bin\tLower Bin" %(sigma)
	for hist in histlist:
		hist = DataFile.Get(hist)
		hist.Scale(1/hist.Integral())
		Sigma(hist,sigma)

print "\nDone!"
