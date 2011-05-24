print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

print "Setting Initial Parameters."
can = TCanvas("Plots","Plots",850,600)
HistogramNames=["PixelBarrelLineardZ","PixelFwdLineardZ","TECLineardZ","TIBLineardZ","TIDLineardZ","TOBLineardZ"]
HistogramNames+=["PixelBarrelNewPVdZ","PixelFwdNewPVdZ","TECNewPVdZ","TIBNewPVdZ","TIDNewPVdZ","TOBNewPVdZ"]
filename="UnweightedTest.root"
pwd = "/data/ndpc2/c/HiggsGammaGamma/CMSSW_4_1_4/src/ND_Hto2Photons/TreeReaders/"
file=TFile(pwd+filename)

for j in range(len(HistogramNames)):
	hist = file.Get(HistogramNames[j])
	print "Looking at Histogram %s" %HistogramNames[j]
	print "%s Total: %i" %(HistogramNames[j],int(hist.GetEntries()))
	print "%s 10mm: %i" %(HistogramNames[j],int(hist.Integral(41,60)))
	print "%s 5mm: %i" %(HistogramNames[j],int(hist.Integral(46,55)))
	print "%s 3mm: %i" %(HistogramNames[j],int(hist.Integral(48,53)))

print "Done"
