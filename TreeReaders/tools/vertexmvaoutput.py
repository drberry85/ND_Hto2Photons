print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)
gStyle.SetCanvasBorderMode(0);
gStyle.SetCanvasColor(kWhite);

print "Setting Initial Parameters."
can = TCanvas("Plots","Plots",750,600)
colors = [kRed+2, kGreen+2]
leg = TLegend(0.50, 0.55, 0.78, 0.85)
leg.SetFillColor(0)
leg.SetBorderSize(1)
leg.SetLineColor(kWhite)
mytext = TLatex()
mytext.SetTextSize(0.04)
mytext.SetNDC()
Legends = ["Non-Jet Tagged Vertex","Jet Tagged Vertex"]
HistNames = [["MVAValuebad","MVAValue"],["MVAResbad","MVARes"]]
FillStyle = [3004,3005]
pwd = "../"
MC = TFile(pwd + "Vertex_Background.root")
Data = TFile(pwd + "Vertex_Data.root")

for i in range(len(HistNames)):
	ScaleFactor=0
	for j in range(len(HistNames[i])):
		DataHist = Data.Get(HistNames[i][j])
		DataHist.Scale(1/DataHist.Integral())
		DataHist.SetLineColor(colors[j])
		DataHist.SetMarkerStyle(20)
		DataHist.SetMarkerColor(colors[j])
		DataHist.SetMaximum(1.1);
		MCHist = MC.Get(HistNames[i][j])
		if (HistNames[i][j]=="MVAResbad"): ScaleFactor=MCHist.Integral()
		MCHist.Scale(1/MCHist.Integral())
		MCHist.SetLineColor(colors[j])
		MCHist.SetFillColor(colors[j])
		MCHist.SetFillStyle(FillStyle[j])
		MCHist.SetMaximum(1.1)
		if (i==0 and j==0):
			DataHist.SetTitle(";BDT Output;")
			DataHist.Draw("e")
		if (i==1 and j==0):
			DummyHist=TH1F("dummy",";MVA_{event};",1,-1.1,1)
			DummyHist.SetMaximum(1.1)
			DummyHist.Draw("AXIS")
			DataHist.Draw("esame")
		else: DataHist.Draw("esame")
		MCHist.Draw("hist same")
		leg.AddEntry(DataHist,Legends[j]+" DATA")
		leg.AddEntry(MCHist,Legends[j]+" MC")
		if (HistNames[i][j]=="MVAResbad"):
			SimHist = MC.Get("MVAResgoodbadsim")
			SimHist.Scale(1/ScaleFactor)
			SimHist.SetLineColor(kRed+3)
			SimHist.SetLineWidth(2)
			SimHist.Draw("hist same")
			leg.AddEntry(SimHist,"Mis-Tagged Vertex MC")
	leg.Draw()
	mytext.DrawLatex(0.47,0.85,"#sqrt{s} = 8 TeV Run 2012A+B+C+D")
	can.SetLogy(0)
	can.SaveAs(HistNames[i][j]+".png")
	can.SetLogy(1)
	can.SaveAs(HistNames[i][j]+"_log.png")
	can.Clear()
	leg.Clear()
	gDirectory.Delete("*")

print "Done!"
