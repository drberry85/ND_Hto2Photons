print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
import math
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

print "Setting Initial Parameters."
can = TCanvas("Plots","Plots",850,650)
can.SetGrid()
leg = TLegend(0.7, 0.7, 0.9, 0.9)
leg.SetFillColor(0)
leg.SetBorderSize(1)
HistogramName="ZMass"
Categories=["","_cat0","_cat1","_cat2","_cat3"]
pwd = "/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/"
#File = TFile(pwd+"ZMuMu_Data_NoCuts.root")
File = TFile(pwd+"ZMuMu_Fall11_NoCuts.root")

for Cat in Categories:
	TotalHist=File.Get(HistogramName+Cat)
	TotalErr=Double(0.0)
	Total=TotalHist.IntegralAndError(1,TotalHist.GetNbinsX(),TotalErr)
	PassHist=File.Get("Pass"+HistogramName+Cat)
	PassErr=Double(0.0)
	Pass=PassHist.IntegralAndError(1,PassHist.GetNbinsX(),PassErr)
	FailHist=File.Get("Fail"+HistogramName+Cat)
	FailErr=Double(0.0)
	Fail=FailHist.IntegralAndError(1,FailHist.GetNbinsX(),FailErr)
	Eff=Pass/Total
	EffErr=Eff*math.sqrt(pow(TotalErr/Total,2)+pow(PassErr/Pass,2))
	print "%s Pass: %0.1f Fail: %0.1f" %(HistogramName+Cat,Pass,Fail)
	print "Efficiency for %s: %0.3f +/- %0.3f" %(HistogramName+Cat,Eff,EffErr)

print "Done!"
