print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

print "Setting Initial Parameters."
can = TCanvas("Plots","Plots",850,650)
colors = [3, 2, 5, 7, 6]
leg = TLegend(0.7, 0.7, 0.9, 0.9)
leg.SetFillColor(0)
leg.SetBorderSize(1)
HistogramNames=["ZMass","ZMassZoom","Numvtx","PhotonEt","PhotonEta","PhotonPhi"]
Categories=["","_cat0","_cat1","_cat2","_cat3"]
FileNames=["ZMuMu_Fall11_NoCuts.root","ZMuMu_TTJets_NoCuts.root"]
Legends = ["Z#mu#mu#gamma","TT+Jets"]
pwd = "/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/"
DataFile = TFile(pwd+"ZMuMu_Data_NoCuts.root")

MCFiles=[]
for file in FileNames:
	MCFiles.append(TFile(pwd+file))

for histbase in HistogramNames:
	for cut in ["","Pass","Fail"]:
		for cat in Categories:
			hist=cut+histbase+cat
			Data=DataFile.Get(hist)
			stack = THStack(hist,"")
			mchistlist=[]
			for i in range(len(MCFiles)):
				MCFiles[i].cd()
				mchistlist.append(MCFiles[i].Get(hist))
				mchistlist[i].SetTitle("")
				if i==0: stack.SetTitle(histbase+";"+mchistlist[0].GetXaxis().GetTitle()+";")
				mchistlist[i].SetLineWidth(0)
				mchistlist[i].SetLineColor(colors[i])
				mchistlist[i].SetFillColor(colors[i])
				stack.Add(mchistlist[i])
				leg.AddEntry(mchistlist[i],Legends[i])
			stack.Draw('hist')
			Data.SetMarkerStyle(20)
			if stack.GetMaximum()<(Data.GetMaximum()+sqrt(Data.GetMaximum())): stack.SetMaximum(Data.GetMaximum()+sqrt(Data.GetMaximum()))
			else: stack.SetMaximum(stack.GetMaximum())
			leg.AddEntry(Data,"Data")
			Data.Draw("esame")
			leg.Draw()
			can.SaveAs("ZMuMuStack_"+hist+"_NoCuts.png")
			stack.Clear()
			leg.Clear()
			gDirectory.Delete("*")
