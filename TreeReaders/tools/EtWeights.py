#print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
import array
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

DataFile=TFile("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/Vertex_Data.root")
MCFile=TFile("/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/Vertex_PJets.root")
DataHist=DataFile.Get("PhotonPtWeight")
MCHist=MCFile.Get("PhotonPtWeight")
#DataHist=DataFile.Get("Numvtx")
#MCHist=MCFile.Get("Numvtx")
DataHist.Scale(1/DataHist.Integral())
MCHist.Scale(1/MCHist.Integral())

for i in range(1,DataHist.GetNbinsX()+2):
	if DataHist.GetBinContent(i)==0 or MCHist.GetBinContent(i)==0:
		print "EtMap[%i]=%f;" %(i-1,0)
	else:
		print "EtMap[%i]=%f;" %(i-1,DataHist.GetBinContent(i)/MCHist.GetBinContent(i))

#print "Done"
