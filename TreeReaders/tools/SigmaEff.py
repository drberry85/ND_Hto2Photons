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
filename="Vertex_Data.root"
pwd = "/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/"
file=TFile(pwd+filename)

for histlist in HistogramNames:
	names=""
	results=""
	resultsRMS=""
	for histname in histlist:
		hist=file.Get(histname)
		RMS=hist.GetRMS()
		sigma = effSigma(hist)
		results+=histname+": "+str(round(sigma,3))[:5]+" "
		resultsRMS+=histname+": "+str(round(RMS,3))[:5]+" "
		#print "%s: %f" %(histname,sigma)
	print "Effective Sigma: "+results
	#print "RMS: \t\t"+resultsRMS

print "Done!"
