print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
import os
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

can = TCanvas("Plots","Plots",850,650)
leg = TLegend(0.6, 0.6, 0.9, 0.9)
leg.SetFillColor(0)
leg.SetBorderSize(1)
pwd="../"
FileList=["Vertex_Data.root","Vertex_Background.root","Vertex_DiPhotonJets.root","Vertex_PJets.root","Vertex_QCD.root","Vertex_WJets.root","Vertex_ZJets.root"]
colors = [1, 2, 4, 7, 3, 6, 5, 9, 8]
file=[]

for i in range(len(FileList)):
	print "Looking at file: %s" %(FileList[i])
	file.append(TFile(pwd+FileList[i]))
	hist=file[i].Get("Numvtx")
	hist.GetXaxis().SetRangeUser(0,50)
	hist.Scale(1/hist.Integral())
	hist.SetLineWidth(2)
	hist.SetLineColor(colors[i])
	hist.SetMaximum(0.1)
	if (i==0): hist.Draw()
	else: hist.Draw("same")
	leg.AddEntry(hist,file[i].GetName()[file[i].GetName().rfind("/")+1:])

leg.Draw()
can.SaveAs("VertexPileup.pdf")
print "Done!"
