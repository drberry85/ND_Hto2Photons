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
colors=[2, 1, 3, 7, 4, 6, 5, 9, 8]
HistogramNames=["ConvR","ConvR-0.5","ConvR-0.8","ConvR-1.0","ConvR-1.3","ConvR-1.5","ConvR0.0","ConvR0.5","ConvR0.8","ConvR1.0","ConvR1.3","ConvRFine-0.5","ConvRFine-0.8","ConvRFine-1.0","ConvRFine-1.3","ConvRFine-1.5","ConvRFine0.0","ConvRFine0.5","ConvRFine0.8","ConvRFine1.0","ConvRFine1.3"]
FileNames=["Vertex_Data.root","Vertex_PhotonPlusJetMC.root"]
leg = TLegend(0.7, 0.7, 0.9, 0.9)
leg.SetFillColor(0)
leg.SetBorderSize(1)
Legends = ["Run2011A+B","#gamma+Jet MC"]
pwd = "/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/"

files=[]
for filename in FileNames:
	files.append(TFile(pwd+filename))

for histname in HistogramNames:
	hist=[]
	for i in range(len(files)):
		files[i].cd()
		hist.append(files[i].Get(histname))
		hist[i].Sumw2()
		hist[i].Scale(1/hist[i].Integral())
		hist[i].SetLineColor(colors[i])
		hist[i].SetLineWidth(2)
		leg.AddEntry(hist[i],Legends[i])
		hist[i].GetXaxis().SetTitle("Radius (cm)")
		hist[i].GetYaxis().SetTitle("Event Fraction")
		hist[i].SetMaximum(0.04)
		if i==0: hist[i].Draw("hist e")
		else: hist[i].Draw("hist e same")
	leg.Draw()
	can.SaveAs(histname+".png")
	can.Clear()
	if histname=="ConvR":
		RatioHist=hist[0]
		RatioHist.Divide(hist[1])
		RatioHist.SetNameTitle("RatioHist",";Conversion Radius (cm);Data/MC Ratio");
		RatioHist.SetMarkerStyle(20)
		RatioHist.SetMarkerSize(1)
		RatioHist.SetLineColor(kBlack)
		RatioHist.SetMinimum(0)
		RatioHist.SetMaximum(4.5)
		RatioHist.Draw("")
		can.SaveAs(histname+"_Ratio.png")
		can.Clear()
	leg.Clear()
	gDirectory.Delete("*")

print "Done"
