print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)
gStyle.SetCanvasBorderMode(0);
gStyle.SetCanvasColor(kWhite);


print "Setting Initial Parameters."
can = TCanvas("Plots","Plots",850,600)
#can.SetLogy(1)
colors = [4, 2, 5, 7, 3, 6]
HistogramNames=["goodconvchi2Barrel","goodconvchi2Endcap","goodconvdistminBarrel","goodconvdistminEndcap"]
HistogramNames+=["goodconveopBarrel","goodconveopEndcap","goodconvrBarrel","goodconvrEndcap"]
HistogramNames+=["goodconvdetaBarrel","goodconvdetaEndcap","goodconvdphiBarrel","goodconvdphiEndcap"]
HistogramNames+=["CosThetaStar_2gammaGoldenSelBarrel","CosThetaStar_2gammaGoldenSelEndcap","pt_2gammaGoldenSelBarrel","pt_2gammaGoldenSelEndcap"]
HistogramNames+=["mass_2gammaSelBarrel","mass_2gammaSelEndcap"]
HistogramNames+=["mass_2gammaGoldenSelBarrel","mass_2gammaGoldenSelEndcap"]
HistogramNames+=["mass_2gammaleftoverSelBarrel","mass_2gammaleftoverSelEndcap"]
Filenames=["QCDDoubleEMEnriched.root","PhotonJetEMEnriched.root","Born.root","Box.root","DrellYan.root"]
file=[]
leg = TLegend(0.6, 0.7, 0.9, 0.9)
leg.SetFillColor(0)
leg.SetBorderSize(1)
Legends = ["QCD","Photon Plus Jet","Born","Box","DrellYan"]
pwd = "/data/ndpc2/c/HiggsGammaGamma/CMSSW_4_1_4/src/ND_Hto2Photons/TreeReaders/"
#Data = TFile(pwd + "Data.root")
#lumi = 94.4
Data = TFile(pwd + "PromptReco.root")
lumi = 228.6

for i in range(len(Filenames)):
	file.append(TFile(pwd+Filenames[i]))

for i in range(len(HistogramNames)):
	DataHist = Data.Get(HistogramNames[i])
	hist=[]
	stack = THStack(HistogramNames[i],"")
	for j in range(len(Filenames)):
		file[j].cd()
		hist.append(file[j].Get(HistogramNames[i]))
		hist[j].Scale(lumi/1000)
		hist[j].SetTitle("");
		if (j==0): stack.SetTitle(hist[0].GetTitle()+";"+hist[0].GetXaxis().GetTitle()+";")
		hist[j].SetLineWidth(0)
		hist[j].SetLineColor(colors[j])
		hist[j].SetFillColor(colors[j])
		stack.Add(hist[j])
		leg.AddEntry(hist[j],Legends[j])
	stack.SetMinimum(0.1)
	stack.Draw()
	if (DataHist.GetMaximum()>stack.GetMaximum()):
		stack.SetMaximum(DataHist.GetMaximum()+sqrt(DataHist.GetMaximum()))
        DataHist.SetMarkerStyle(20);
	leg.AddEntry(DataHist,"Data")
	DataHist.Draw("esame")
	leg.Draw()
#	can.SaveAs(HistogramNames[i]+"DataLogScale.png")
	can.SaveAs(HistogramNames[i]+"Data.png")
	stack.Clear()
	leg.Clear()
	gDirectory.Delete("*")
