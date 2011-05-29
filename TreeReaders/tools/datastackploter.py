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
HistogramNames+=["goodconvzBarrel","goodconvzEndcap","goodconveta"]
HistogramNames+=["CosThetaStar_2gammaSelBarrel","CosThetaStar_2gammaSelEndcap"]
HistogramNames+=["CosThetaStar_2gammaGoldenSelBarrel","CosThetaStar_2gammaGoldenSelEndcap"]
HistogramNames+=["CosThetaStar_2gammaleftoverSelBarrel","CosThetaStar_2gammaleftoverSelEndcap"]
HistogramNames+=["mass_2gammaSelBarrel","mass_2gammaSelEndcap"]
HistogramNames+=["mass_2gammaGoldenSelBarrel","mass_2gammaGoldenSelEndcap"]
HistogramNames+=["mass_2gammaleftoverSelBarrel","mass_2gammaleftoverSelEndcap"]
Normalization = [0.0]*len(HistogramNames)
Filenames=["QCDDoubleEMEnriched.root","PhotonJetEMEnriched.root","Born.root","Box.root","DrellYan.root"]
file=[]
leg = TLegend(0.7, 0.7, 0.9, 0.9)
leg.SetFillColor(0)
leg.SetBorderSize(1)
Legends = ["QCD","Photon Plus Jet","Born","Box","DrellYan"]
pwd = "/data/ndpc2/c/HiggsGammaGamma/CMSSW_4_1_4/src/ND_Hto2Photons/TreeReaders/"
#Data = TFile(pwd + "Data.root")
#lumi = 94.4
Data = TFile(pwd + "PromptReco.root")
lumi = 229.8

for i in range(len(Filenames)):
	file.append(TFile(pwd+Filenames[i]))
	for j in range(len(HistogramNames)):
		hist=file[i].Get(HistogramNames[j])
		Normalization[j] += hist.Integral()

for i in range(len(HistogramNames)):
	DataHist = Data.Get(HistogramNames[i])
	DataNorm = DataHist.Integral()
	hist=[]
	stack = THStack(HistogramNames[i],"")
	for j in range(len(Filenames)):
		file[j].cd()
		hist.append(file[j].Get(HistogramNames[i]))
		#hist[j].Scale(lumi/1000)
		hist[j].Scale(DataNorm/Normalization[i])
		hist[j].SetTitle("");
		if (j==0): stack.SetTitle(hist[0].GetTitle()+";"+hist[0].GetXaxis().GetTitle()+";")
		hist[j].SetLineWidth(0)
		hist[j].SetLineColor(colors[j])
		hist[j].SetFillColor(colors[j])
		stack.Add(hist[j])
		leg.AddEntry(hist[j],Legends[j])
	if (can.GetLogy()==1): stack.SetMinimum(0.1)
	stack.Draw()
	if (HistogramNames[i].find("mass")!=-1): stack.GetXaxis().SetRangeUser(100,160)
	DataHist.SetMarkerStyle(20);
	if (stack.GetMaximum()<(DataHist.GetMaximum()+sqrt(DataHist.GetMaximum()))): stack.SetMaximum(DataHist.GetMaximum()+sqrt(DataHist.GetMaximum()))
	else: stack.SetMaximum(stack.GetMaximum())
	if can.GetLogy()==0 and (HistogramNames[i].find("CosThetaStar")!=-1 or HistogramNames[i].find("goodconveta")!=-1): stack.SetMaximum(1.5*(DataHist.GetMaximum()+sqrt(DataHist.GetMaximum())))
	if can.GetLogy()==1 and (HistogramNames[i].find("CosThetaStar")!=-1 or HistogramNames[i].find("mass")!=-1 or HistogramNames[i].find("goodconveta")): stack.SetMaximum(10*(DataHist.GetMaximum()+sqrt(DataHist.GetMaximum())))
	leg.AddEntry(DataHist,"Data")
	DataHist.Draw("esame")
	leg.Draw()
	if (can.GetLogy()==1): outputname=HistogramNames[i]+"PromptRecoLogScale.png"
	else: outputname=HistogramNames[i]+"PromptReco.png"
	can.SaveAs(outputname)
	stack.Clear()
	leg.Clear()
	gDirectory.Delete("*")
