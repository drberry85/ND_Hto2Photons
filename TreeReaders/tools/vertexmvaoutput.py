print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)
gStyle.SetCanvasBorderMode(0);
gStyle.SetCanvasColor(kWhite);

print "Setting Initial Parameters."
can = TCanvas("Plots","Plots",850,600)
can.SetLogy()
colors = [kRed+2, kGreen+2]
leg = TLegend(0.65, 0.55, 0.85, 0.75)
leg.SetFillColor(0)
leg.SetBorderSize(1)
leg.SetLineColor(kWhite)
mytext = TLatex()
mytext.SetTextSize(0.04)
mytext.SetNDC()
Legends = ["Wrong Vertex","Right Vertex"]
HistNames = ["MVAValuebad","MVAValue"]
FillStyle = [3004,3005]
pwd = "/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/"
MC = TFile(pwd + "Vertex_PhotonPlusJetMC.root")
Data = TFile(pwd + "Vertex_Data.root")

for i in range(len(HistNames)):
	DataHist = Data.Get(HistNames[i])
	DataHist.Sumw2()
	DataHist.Scale(1/DataHist.Integral())
	DataHist.SetLineColor(colors[i])
	DataHist.SetMarkerStyle(20)
	DataHist.SetMarkerColor(colors[i])
	MCHist = MC.Get(HistNames[i])
	MCHist.Scale(1/MCHist.Integral())
	MCHist.SetLineColor(colors[i])
	MCHist.SetFillColor(colors[i])
	MCHist.SetFillStyle(FillStyle[i])
	if (i==0):
		DataHist.SetTitle(";BDT Output;")
		DataHist.Draw("e")
	else: DataHist.Draw("esame")
	MCHist.Draw("same")
	leg.AddEntry(DataHist,Legends[i]+" DATA")
	leg.AddEntry(MCHist,Legends[i]+" MC")

leg.Draw()
mytext.DrawLatex(0.55,0.8,"#sqrt{s} = 7 TeV Run 2011A+B")
outputname="MVAOutput_log.png"
can.SaveAs(outputname)
