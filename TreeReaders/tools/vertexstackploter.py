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
colors = [3, 2, 5, 7, 6]
HistogramNames=["PhotonPt","ConvPt","PixelBarrelConvdZ","PixelFwdConvdZ","TECConvdZ","TIBConvdZ","TIDConvdZ","TOBConvdZ"]
HistogramNames+=["PixelBarrelSuperdZ","PixelFwdSuperdZ","TECSuperdZ","TIBSuperdZ","TIDSuperdZ","TOBSuperdZ"]
Normalization = [0.0]*len(HistogramNames)
Filenames=["Vertex_PhotonPlusJetMC.root","Vertex_QCDEMEnriched.root"]
file=[]
leg = TLegend(0.7, 0.7, 0.9, 0.9)
leg.SetFillColor(0)
leg.SetBorderSize(1)
Legends = ["#gamma+Jet","QCD EMEnriched"]
pwd = "/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/"
Data = TFile(pwd + "Vertex_Data.root")

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
		hist[j].Scale(DataNorm/Normalization[i])
		hist[j].SetTitle("");
		if j==0: stack.SetTitle(hist[0].GetTitle()+";"+hist[0].GetXaxis().GetTitle()+";")
		hist[j].SetLineWidth(0)
		hist[j].SetLineColor(colors[j])
		hist[j].SetFillColor(colors[j])
		stack.Add(hist[j])
		leg.AddEntry(hist[j],Legends[j])
	if can.GetLogy()==1: stack.SetMinimum(0.1)
	stack.Draw()
	if HistogramNames[i].find("mass")!=-1: stack.GetXaxis().SetRangeUser(100,160)
	DataHist.SetMarkerStyle(20)
	if stack.GetMaximum()<(DataHist.GetMaximum()+sqrt(DataHist.GetMaximum())): stack.SetMaximum(DataHist.GetMaximum()+sqrt(DataHist.GetMaximum()))
	else: stack.SetMaximum(stack.GetMaximum())
	leg.AddEntry(DataHist,"Data")
	DataHist.Draw("esame")
	leg.Draw()
	if can.GetLogy()==1: outputname=HistogramNames[i]+"VertexLogScale.png"
	else: outputname=HistogramNames[i]+"Vertex.png"
	can.SaveAs(outputname)
	stack.Clear()
	leg.Clear()
	gDirectory.Delete("*")