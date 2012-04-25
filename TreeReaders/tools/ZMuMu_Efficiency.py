print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

print "Setting Initial Parameters."
can = TCanvas("Plots","Plots",850,650)
can.SetGrid()
leg = TLegend(0.7, 0.7, 0.9, 0.9)
leg.SetFillColor(0)
leg.SetBorderSize(1)
HistogramNames=["Numvtx","PhotonEt","PhotonEta","PhotonPhi"]
Categories=["","_cat0","_cat1","_cat2","_cat3"]
pwd = "/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/"
DataFile = TFile(pwd+"ZMuMu_Data.root")
MCFile = TFile(pwd+"ZMuMu_Fall11.root")

for histbase in HistogramNames:
	for cat in Categories:
		hist=histbase+cat
		Data = DataFile.Get(hist)
		DataPass = DataFile.Get("Pass"+hist)
		MC = MCFile.Get(hist)
		MCPass = MCFile.Get("Pass"+hist)
		DataEff = DataPass.Clone()
		DataEff.Divide(DataPass,Data,1,1,"B")
		DataEff.SetMaximum(1.01)
		DataEff.SetMinimum(0.80)
		MCEff = MCPass.Clone()
		MCEff.Divide(MCPass,MC,1,1,"B")
		DataEff.Draw("")
		if hist.find("Numvtx")!=-1: DataEff.GetXaxis().SetRangeUser(0,40)
		MCEff.SetLineColor(kRed)
		MCEff.Draw("same")
		can.SaveAs("ZMuMuEff_"+hist+".png")
