print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
import array
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

print "Setting Initial Parameters."
can = TCanvas("Plots","Plots",850,650)
can.SetGrid()
leg = TLegend(0.4, 0.9, 0.7, 0.8)
leg.SetFillColor(0)
leg.SetBorderSize(1)
NumberatorHists=["ConvEta","ConvPt"]
DenominatorHists=["PhotonEta","PhotonPt"]
Legends=["Data RunA+B+C+D","PJet MC"]
OutFileNames=["ConversionEtaEff.png","ConversionPtEff.png"]
Maximum=[80,80]
FileNames=["Vertex_Data.root","Vertex_Background.root"]
pwd = "../"
colors=[3, 2, 1, 7, 4, 6, 5, 9, 8]

files=[]
for filename in FileNames:
	files.append(TFile(pwd+filename))

for numhist,denhist,outfile,max in zip(NumberatorHists,DenominatorHists,OutFileNames,Maximum):
	top=[]
	bottom=[]
	for i in range(len(files)):
		top.append(files[i].Get(numhist))
		bottom.append(files[i].Get(denhist))
		top[i].Sumw2()
		bottom[i].Sumw2()
		top[i].Divide(bottom[i])
		top[i].Scale(100)
		top[i].SetLineColor(colors[i])
		top[i].SetLineWidth(2)
		top[i].GetYaxis().SetTitle("Efficiency")
		leg.AddEntry(top[i],Legends[i])
		if i==0:
			top[i].SetMinimum(0)
			top[i].SetMaximum(max)
			top[i].Draw("e")
		else:
			top[i].Draw("esame")
	leg.Draw()
	can.SaveAs(outfile)
	can.Clear()
	if (len(top)==2):
		RatioHist=top[0]
		RatioHist.Divide(top[1])
		RatioHist.SetNameTitle("RatioHist",";"+top[0].GetXaxis().GetTitle()+";Data/MC Ratio");
		RatioHist.SetMarkerStyle(20)
		RatioHist.SetMarkerSize(1)
		RatioHist.SetLineColor(kBlack)
		if numhist.find("Eta") != -1:
			RatioHist.SetMinimum(0.0)
			RatioHist.SetMaximum(2.0)
			RatioHist.SetNameTitle("RatioHist",";#eta;Data/MC Ratio");
		else:
			RatioHist.SetMinimum(0.0)
			RatioHist.SetMaximum(2.0)
			RatioHist.SetNameTitle("RatioHist",";Pt (GeV);Data/MC Ratio");
		RatioHist.Draw("")
		can.SaveAs(outfile.replace(".png","_Ratio.png"))
		can.Clear()
	leg.Clear()
	gDirectory.Delete("*")

print "Done"
