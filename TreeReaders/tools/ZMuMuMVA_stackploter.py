print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
#gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

print "Setting Initial Parameters."
can = TCanvas("Plots","Plots",850,650)
can.SetLogy()
colors = [3, 2, 5, 7, 6]
leg = TLegend(0.7, 0.8, 0.9, 0.9)
leg.SetFillColor(0)
leg.SetBorderSize(1)
HistogramNames=["pho_pt","pho_idmvanew","pho_tmva_photonid_etawidth","pho_tmva_photonid_phiwidth","pho_tmva_photonid_r9","pho_tmva_photonid_sceta","pho_tmva_photonid_sieie","pho_tmva_photonid_sieip","pho_tmva_photonid_s4ratio","pho_tmva_photonid_lambdaratio","pho_tmva_photonid_eventrho","pho_tmva_photonid_ESEffSigmaRR","pho_tmva_photonid_pfchargedisogood03","pho_tmva_photonid_pfchargedisobad03","pho_tmva_photonid_pfphotoniso03","pho_tmva_photonid_pfneutraliso03"]
Categories=["Barrel","Endcap"]
pwd = "/data/pccmsnd1/b/ZMuMuGam/CMSSW_4_4_0/src/ND_Hto2Photons/TreeReaders/2012RootFiles/new2012GB_powheg/"
pwddata = "/data/pccmsnd1/b/ZMuMuGam/CMSSW_4_4_0/src/ND_Hto2Photons/TreeReaders/2012RootFiles/new2012GB_powheg/"
FileNames=["ZMuMu_Summer12_PhoPFPresel_HighPt.root"]
#FileNames=["ZMuMu_Summer12_HighPt.root"]
Legends = ["Z#mu#mu#gamma"]
DataFile = TFile(pwddata+"ZMuMu_Run2012_PhoPFPresel_HighPt.root")
#DataFile = TFile(pwddata+"ZMuMu_Run2012_HighPt.root")


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
		stack.Draw('histe')
		Data.SetMarkerStyle(20)
		if stack.GetMaximum()<(Data.GetMaximum()+sqrt(Data.GetMaximum())): stack.SetMaximum(Data.GetMaximum()+sqrt(Data.GetMaximum()))
		else: stack.SetMaximum(stack.GetMaximum())
                stack.SetMinimum(0.1)
		leg.AddEntry(Data,"Data")
		Data.Draw("esame")
		leg.Draw()
		can.SaveAs("ZMuMuMVAStack_"+hist+".png")
		RatioHist = Data.Clone("")
		RatioHist.Divide(Data,MCHistSum,1,1,"B")
		RatioHist.SetMinimum(-2)
		RatioHist.SetMaximum(10)
		can.SetLogy(0)
		RatioHist.Draw("")
		MyLine = TLine(RatioHist.GetXaxis().GetXmin(),1,RatioHist.GetXaxis().GetXmax(),1)
		MyLine.Draw("")
		can.SaveAs("ZMuMuMVAStack_"+hist+"_Ratio.png")
		can.SetLogy()
		stack.Clear()
		leg.Clear()
		gDirectory.Delete("*")
