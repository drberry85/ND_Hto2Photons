print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

print "Setting Initial Parameters."
can = TCanvas("Plots","Plots",850,600)
HistogramNames=["convphotonr9Barrel","convphotonr9Endcap","goodconvphotonr9Barrel","goodconvphotonr9Endcap"]
filename="Conversion_UnweightedHiggsAnalysis115GeV.root"
pwd = "/data/ndpc2/c/HiggsGammaGamma/CMSSW_4_1_4/src/ND_Hto2Photons/TreeReaders/"
file=TFile(pwd+filename)

for j in range(len(HistogramNames)):
	hist = file.Get(HistogramNames[j])
	print "Looking at Histogram %s" %HistogramNames[j]
	print "%s R9>0.93 %i" %(HistogramNames[j],int(hist.Integral(94,100)))
	print "%s R9<0.93 %i" %(HistogramNames[j],int(hist.Integral(1,93)))

print "Done"
