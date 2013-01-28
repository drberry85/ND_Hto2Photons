print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

print "Setting Initial Parameters."
can = TCanvas("Plots","Plots",850,650)
colors = [3, 2, 5, 7, 6]
leg = TLegend(0.7, 0.7, 0.9, 0.9)
leg.SetFillColor(0)
leg.SetBorderSize(1)
HistogramNames=["ZMass","ZMassZoom","Numvtx","PhotonEt","PhotonEta","PhotonPhi"]
Categories=["","_cat0","_cat1","_cat2","_cat3"]
FileNames=["ZMuMu_ZToMuMu_PhoPFPresel_Corr2012_HighPt.root"]
Legends = ["Simulation"]
pwd = "../"
DataFile = TFile(pwd+"ZMuMu_Data_PhoPFPresel_HighPt.root")

MCFiles=[]
for file in FileNames:
    MCFiles.append(TFile(pwd+file))

for histbase in HistogramNames:
    for cut in ["","Pass","Fail"]:
        for cat in Categories:
            hist=cut+histbase+cat
            Data=DataFile.Get(hist)
            stack = THStack(hist,"")
            mchistlist=[]
            Normalization=0
            for i in range(len(MCFiles)): Normalization+=MCFiles[i].Get(hist).Integral()
            for i in range(len(MCFiles)):
                MCFiles[i].cd()
                mchistlist.append(MCFiles[i].Get(hist))
                mchistlist[i].Scale(Data.Integral()/Normalization)
                mchistlist[i].SetTitle("")
                if i==0: stack.SetTitle(histbase+";"+mchistlist[0].GetXaxis().GetTitle()+";")
                mchistlist[i].SetLineWidth(0)
                mchistlist[i].SetLineColor(colors[i])
                mchistlist[i].SetFillColor(colors[i])
                stack.Add(mchistlist[i])
                leg.AddEntry(mchistlist[i],Legends[i])
            stack.SetTitle(";"+mchistlist[0].GetXaxis().GetTitle()+";Events in Data")
            stack.Draw('hist')
            Data.SetMarkerStyle(20)
            if stack.GetMaximum()<(Data.GetMaximum()+sqrt(Data.GetMaximum())): stack.SetMaximum(Data.GetMaximum()+sqrt(Data.GetMaximum()))
            else: stack.SetMaximum(stack.GetMaximum())
            leg.AddEntry(Data,"8TeV Data, 19.0 fb^{-1}")
            Data.Draw("esame")
            leg.Draw()
            can.SaveAs("ZMuMuStack_"+hist+"_PhoPFPresel_HighPt.png")
            stack.Clear()
            leg.Clear()
            gDirectory.Delete("*")
