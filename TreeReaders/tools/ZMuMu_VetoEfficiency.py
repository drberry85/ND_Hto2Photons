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
pwd = "/data/pccmsnd1/b/ZMuMuGam/CMSSW_4_4_0/src/ND_Hto2Photons/TreeReaders/2012RootFiles/new2012GB_powheg/allDataJun8/electronveto/"
dataFile = TFile(pwd+"ZMuMu_Run2012_PhoPFPresel_HighPt.root")
mcFile = TFile(pwd+"ZMuMu_Summer12powheg_PhoPFPresel_Corr2012_HighPt.root")

#dataFile = TFile(pwd+"ZMuMu_Run2012_HighPt.root")
#mcFile = TFile(pwd+"ZMuMu_Summer12powheg_Corr2012_HighPt.root")


print " Efficiency in data "
for Cat in Categories:
	TotalHist=dataFile.Get(HistogramName+Cat)
	TotalErr=Double(0.0)
	Total=TotalHist.IntegralAndError(1,TotalHist.GetNbinsX(),TotalErr)
	PassHist=dataFile.Get("Pass"+HistogramName+Cat)
	PassErr=Double(0.0)
	Pass=PassHist.IntegralAndError(1,PassHist.GetNbinsX(),PassErr)
	FailHist=dataFile.Get("Fail"+HistogramName+Cat)
	FailErr=Double(0.0)
	Fail=FailHist.IntegralAndError(1,FailHist.GetNbinsX(),FailErr)
	Eff=Pass/Total
        EffErr= sqrt(Eff*(1.-Eff)/Total) 
	print "%s Pass: %0.1f Fail: %0.1f  Efficiency for %s: %0.3f +/- %0.3f" %(HistogramName+Cat,Pass,Fail,HistogramName+Cat,Eff,EffErr)
	#print "Efficiency for %s: %0.3f +/- %0.3f" %(HistogramName+Cat,Eff,EffErr)

print " Efficiency in MC"
for Cat in Categories:
	TotalHist=mcFile.Get(HistogramName+Cat)
	TotalErr=Double(0.0)
	Total=TotalHist.IntegralAndError(1,TotalHist.GetNbinsX(),TotalErr)
	PassHist=mcFile.Get("Pass"+HistogramName+Cat)
	PassErr=Double(0.0)
	Pass=PassHist.IntegralAndError(1,PassHist.GetNbinsX(),PassErr)
	FailHist=mcFile.Get("Fail"+HistogramName+Cat)
	FailErr=Double(0.0)
	Fail=FailHist.IntegralAndError(1,FailHist.GetNbinsX(),FailErr)
	Eff=Pass/Total
	EffErr= sqrt(Eff*(1.-Eff)/Total) 
	print "%s Pass: %0.1f Fail: %0.1f  Efficiency for %s: %0.3f +/- %0.3f" %(HistogramName+Cat,Pass,Fail,HistogramName+Cat,Eff,EffErr)
	#print "Efficiency for %s: %0.3f +/- %0.3f" %(HistogramName+Cat,Eff,EffErr)

print "Ratios data/MC"
for Cat in Categories:
	dataTotalHist=dataFile.Get(HistogramName+Cat)
	dataTotalErr=Double(0.0)
	dataTotal=dataTotalHist.IntegralAndError(1,dataTotalHist.GetNbinsX(),dataTotalErr)
	dataPassHist=dataFile.Get("Pass"+HistogramName+Cat)
	dataPassErr=Double(0.0)
	dataPass=dataPassHist.IntegralAndError(1,dataPassHist.GetNbinsX(),dataPassErr)
	dataFailHist=dataFile.Get("Fail"+HistogramName+Cat)
	dataFailErr=Double(0.0)
	dataFail=dataFailHist.IntegralAndError(1,dataFailHist.GetNbinsX(),dataFailErr)
	dataEff=dataPass/dataTotal
	dataEffErr= sqrt(dataEff*(1.-dataEff)/dataTotal) 
	#print "%s Pass: %0.1f Fail: %0.1f  Efficiency for %s: %0.3f +/- %0.3f" %(HistogramName+Cat,dataPass,dataFail,HistogramName+Cat,dataEff,dataEffErr)

	mcTotalHist=mcFile.Get(HistogramName+Cat)
	mcTotalErr=Double(0.0)
	mcTotal=mcTotalHist.IntegralAndError(1,mcTotalHist.GetNbinsX(),mcTotalErr)
	mcPassHist=mcFile.Get("Pass"+HistogramName+Cat)
	mcPassErr=Double(0.0)
	mcPass=mcPassHist.IntegralAndError(1,mcPassHist.GetNbinsX(),mcPassErr)
	mcFailHist=mcFile.Get("Fail"+HistogramName+Cat)
	mcFailErr=Double(0.0)
	mcFail=mcFailHist.IntegralAndError(1,mcFailHist.GetNbinsX(),mcFailErr)
	mcEff=mcPass/mcTotal
	mcEffErr= sqrt(mcEff*(1.-mcEff)/mcTotal) 
	#print "%s Pass: %0.1f Fail: %0.1f  Efficiency for %s: %0.3f +/- %0.3f" %(HistogramName+Cat,mcPass,mcFail,HistogramName+Cat,mcEff,mcEffErr)
        ratio=dataEff/mcEff
	ratioErr=(dataEffErr*dataEffErr)/mcEff*mcEff  + (dataEff*dataEff)/(mcEff*mcEff*mcEff*mcEff) * mcEffErr*mcEffErr
	ratioErr=sqrt(ratioErr)
	print "Ratio %0.3f +- %0.3f" %(ratio,ratioErr)
print "Done!"
