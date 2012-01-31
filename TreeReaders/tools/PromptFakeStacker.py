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
HistogramNames=["PhotonPt","PixelBarrelConvdZ","PixelFwdConvdZ","TECConvdZ","TIBConvdZ","TIDConvdZ","TOBConvdZ"]
HistogramNames+=["PixelBarrelSuperdZ","PixelFwdSuperdZ","TECSuperdZ","TIBSuperdZ","TIDSuperdZ","TOBSuperdZ"]
Normalization = [0.0]*len(HistogramNames)
leg = TLegend(0.7, 0.7, 0.9, 0.9)
leg.SetFillColor(0)
leg.SetBorderSize(1)
pwd = "/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/"
Data = TFile(pwd + "Vertex_PhotonPlusJetData.root")
Prompt = TFile(pwd + "Vertex_Background_Prompt.root")
Fake = TFile(pwd + "Vertex_Background_Fake.root")
lumi = 1078.3

for i in range(len(HistogramNames)):
	DataHist = Data.Get(HistogramNames[i])
	DataNorm = DataHist.Integral()
	PromptHist = Prompt.Get(HistogramNames[i])
	PromptNorm = PromptHist.Integral()
	FakeHist = Fake.Get(HistogramNames[i])
	FakeNorm = FakeHist.Integral()
	stack = THStack(HistogramNames[i],"")

	PromptHist.Scale(DataNorm/(PromptNorm+FakeNorm))
	PromptHist.SetLineColor(9)
	PromptHist.SetFillColor(9)
	stack.SetTitle(PromptHist.GetTitle()+";"+PromptHist.GetXaxis().GetTitle()+";")
	stack.Add(PromptHist)
	leg.AddEntry(PromptHist,"Prompt Photon")

	FakeHist.Scale(DataNorm/(PromptNorm+FakeNorm))
	FakeHist.SetLineColor(8)
	FakeHist.SetFillColor(8)
	stack.Add(FakeHist)
	leg.AddEntry(FakeHist,"Fake Photon")

	stack.Draw()
	DataHist.SetMarkerStyle(20)
	if stack.GetMaximum()<(DataHist.GetMaximum()+sqrt(DataHist.GetMaximum())): stack.SetMaximum(DataHist.GetMaximum()+sqrt(DataHist.GetMaximum()))
	else: stack.SetMaximum(stack.GetMaximum())
	leg.AddEntry(DataHist,"Data")
	DataHist.Draw("esame")
	leg.Draw()
	if can.GetLogy()==1: outputname=HistogramNames[i]+"VertexPromptFakeLogScale.png"
	else: outputname=HistogramNames[i]+"VertexPromptFake.png"
	can.SaveAs(outputname)
	stack.Clear()
	leg.Clear()
	gDirectory.Delete("*")
