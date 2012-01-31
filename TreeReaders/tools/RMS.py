print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
import array
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

print "Setting Initial Parameters."
HistogramNames=[["PixelBarrelSuperdZ","TIBSuperdZ","TOBSuperdZ","PixelFwdSuperdZ","TIDSuperdZ","TECSuperdZ"]]
HistogramNames.append(["PixelBarrelConvdZ","TIBConvdZ","TOBConvdZ","PixelFwdConvdZ","TIDConvdZ","TECConvdZ"])
filename="Vertex_PhotonPlusJetMC.root"
pwd = "/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/"
file=TFile(pwd+filename)

for histlist in HistogramNames:
	names=""
	results=""
	for histname in histlist:
		hist=file.Get(histname)
		RMS = hist.GetRMS()
		results+=histname+": "+str(round(RMS,3))[:5]+" "
	print results

print "Done!"
