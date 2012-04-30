print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

print "Setting Initial Parameters."
can = TCanvas("Plots","Plots",850,650)
can.SetLogy()
colors = [3, 2, 5, 7, 6]
leg = TLegend(0.7, 0.8, 0.9, 0.9)
leg.SetFillColor(0)
leg.SetBorderSize(1)
HistogramNames=["pho_idmva","pho_tmva_id_mit_ecal","pho_tmva_id_mit_etawidth","pho_tmva_id_mit_hcal","pho_tmva_id_mit_hoe","pho_tmva_id_mit_nvtx","pho_tmva_id_mit_phiwidth","pho_tmva_id_mit_preshower","pho_tmva_id_mit_r9","pho_tmva_id_mit_sceta","pho_tmva_id_mit_sieie","pho_tmva_id_mit_tiso1","pho_tmva_id_mit_tiso2","pho_tmva_id_mit_tiso3"]
Categories=["Barrel","Endcap"]
FileNames=["ZMuMu_Fall11.root","ZMuMu_TTJets.root"]
Legends = ["Z#mu#mu#gamma","TT+Jets"]
pwd = "/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/"
DataFile = TFile(pwd+"ZMuMu_Data.root")

MCFiles=[]
for file in FileNames:
	MCFiles.append(TFile(pwd+file))

for histbase in HistogramNames:
	for cat in Categories:
		hist=histbase+cat
		Data=DataFile.Get(hist)
		stack = THStack(hist,"")
		mchistlist=[]
		MCHistSum = TH1F()
		for i in range(len(MCFiles)):
			MCFiles[i].cd()
			mchistlist.append(MCFiles[i].Get(hist))
			if i==0: MCHistSum=mchistlist[i]
			else: MCHistSum.Add(mchistlist[i])
			mchistlist[i].SetTitle("")
			if i==0: stack.SetTitle(hist+";"+mchistlist[0].GetXaxis().GetTitle()+";")
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
		can.SaveAs("ZMuMuMVAStack_"+hist+".png")
		RatioHist = Data.Clone("")
		#import pdb; pdb.set_trace()
		RatioHist.Divide(Data,MCHistSum,1,1,"B")
		can.SetLogy(0)
		RatioHist.Draw("")
		can.SaveAs("ZMuMuMVAStack_"+hist+"_Ratio.png")
		can.SetLogy()
		stack.Clear()
		leg.Clear()
		gDirectory.Delete("*")
