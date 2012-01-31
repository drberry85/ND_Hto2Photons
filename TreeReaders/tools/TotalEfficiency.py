print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
from math import pow
import array
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

print "Setting Initial Parameters."
filenamelist=["Vertex_Data.root","Vertex_PhotonPlusJetMC.root"]
histname="MVAdZ"
pwd = "/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/"
LowerBin=41
UpperBin=60

for filename in filenamelist:
	file = TFile(pwd+filename)
	hist = file.Get(histname)
	NumeratorError = Double(0.0)
	Numerator = hist.IntegralAndError(LowerBin,UpperBin,NumeratorError,"cp")
	DenominatorError = Double(0.0)
	Denominator = hist.IntegralAndError(0,101,DenominatorError,"cp")
	Value = Numerator/Denominator
	Error = Value*sqrt(pow(NumeratorError/Numerator,2)+pow(DenominatorError/Denominator,2))
	print "For File: %s Effeceny is: %0.3f+%0.3f" %(filename,Value,Error)
