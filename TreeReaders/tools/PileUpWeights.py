#print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
import array
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

DataFile=TFile("/data/ndpc2/b/drberry/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/Vertex_Data.root")
DataHist=DataFile.Get("Numvtx")
MCFile=TFile("/data/ndpc2/b/drberry/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/Vertex_PJet_NoPU.root")
MCHist=MCFile.Get("Numvtx")

for i in range(DataHist.GetNbinsX()+2):
	if i==0: continue
	DataHist.Scale(1/DataHist.Integral())
	MCHist.Scale(1/MCHist.Integral())
	if DataHist.GetBinContent(i)==0 or MCHist.GetBinContent(i)==0:
		print "PileUpMap[%i]=%f;" %(i-1,0)
	else:
		print "PileUpMap[%i]=%f;" %(i-1,DataHist.GetBinContent(i)/MCHist.GetBinContent(i))

#print "Done"
