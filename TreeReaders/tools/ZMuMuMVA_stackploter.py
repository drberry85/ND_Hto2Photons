print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
#gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

print "Setting Initial Parameters."
can = TCanvas("Plots","Plots",850,650)
#can.SetLogy()
colors = [3, 2, 5, 7, 6]
datacolors = [1,2,3,4,5,6,7]
leg = TLegend(0.7, 0.8, 0.9, 0.9)
leg.SetFillColor(0)
leg.SetBorderSize(1)
HistogramNames=["pho_pt","pho_mitmva","pho_tmva_photonid_etawidth","pho_tmva_photonid_phiwidth","pho_tmva_photonid_r9","pho_tmva_photonid_sceta","pho_tmva_photonid_sieie","pho_tmva_photonid_sieip","pho_tmva_photonid_s4ratio","pho_tmva_photonid_lambdaratio","pho_tmva_photonid_eventrho","pho_tmva_photonid_ESEffSigmaRR","pho_tmva_photonid_pfchargedisogood03","pho_tmva_photonid_pfchargedisobad03","pho_tmva_photonid_pfphotoniso03","pho_tmva_photonid_pfneutraliso03"]
Categories=["Barrel","Endcap"]
pwd = "../"
DataFileNames = ["ZMuMu_Data_PhoPFPresel_HighPt.root"]
DataLegends = ["8TeV Data, 19.1 fb^{-1}"]
#DataFileNames = ["ZMuMu_Run2012ABC.root","ZMuMu_Run2012D.root"]
#DataLegends = ["Run2012ABC, 11.8fb","Run2012D, 7.3fb scaled to 11.8fb"]
#FileNames = []
FileNames = ["ZMuMu_ZToMuMu_PhoPFPresel_Corr2012_HighPt.root"]
Legends = ["Z#mu#mu#gamma"]
label="_PhoPFPresel_HighPt"
#label="_RunComparisonNoSelection"
MakeLogScale=False
if MakeLogScale: label+="_Log"
MCFiles=[]
DataFiles=[]
for file in FileNames: MCFiles.append(TFile(pwd+file))
for file in DataFileNames: DataFiles.append(TFile(pwd+file))

for histbase in HistogramNames:
    for cat in Categories:
        can.SetCanvasSize(850,650)
        if MakeLogScale: can.SetLogy()
        if histbase=="pho_tmva_photonid_sceta" and cat=="Endcap": hist=histbase
        else: hist=histbase+cat
        stack = THStack(hist,"")
        mchistlist=[]
        datahistlist=[]
        MCHistSum = TH1F()
        Normalization=0
        for i in range(len(MCFiles)): Normalization+=MCFiles[i].Get(hist).Integral()
        for i in range(len(MCFiles)):
            MCFiles[i].cd()
            mchistlist.append(MCFiles[i].Get(hist))
            mchistlist[i].Scale(DataFiles[0].Get(hist).Integral()/Normalization)
            if i==0: MCHistSum=mchistlist[i]
            else: MCHistSum.Add(mchistlist[i])
            mchistlist[i].SetTitle("")
            if i==0: stack.SetTitle(hist+";"+mchistlist[0].GetXaxis().GetTitle()+";")
            mchistlist[i].SetLineWidth(0)
            mchistlist[i].SetLineColor(colors[i])
            mchistlist[i].SetFillColor(colors[i])
            stack.Add(mchistlist[i])
            leg.AddEntry(mchistlist[i],Legends[i])
        if len(MCFiles)>0:
            stack.SetTitle(";"+mchistlist[0].GetXaxis().GetTitle()+";Events in Data")
            stack.SetMaximum(1.3*stack.GetMaximum())
            stack.Draw('hist')
        for i in range(len(DataFiles)):
            datahistlist.append(DataFiles[i].Get(hist))
            datahistlist[i].SetMarkerStyle(20)
            datahistlist[i].SetMarkerColor(datacolors[i])
            leg.AddEntry(datahistlist[i],DataLegends[i])
            if DataLegends[i].find("Run2012D")!=-1: datahistlist[i].Scale(11.8/7.3)
            if i==0: datahistlist[0].SetMaximum(datahistlist[0].GetMaximum()*1.3)
            if i==0 and len(MCFiles)==0: datahistlist[i].Draw("e")
            else: datahistlist[i].Draw("esame")
            leg.Draw()
        can.SaveAs("ZMuMuMVAStack_"+hist+label+".png")
        if len(DataFiles)+len(MCFiles)==2:
            can.SetCanvasSize(850,300)
            can.SetLogy(0)
            if len(DataFiles)==0: RatioHist = mchistlist[0].Clone("")
            else: RatioHist = datahistlist[0].Clone("")
            if len(DataFiles)==1: RatioHist.Divide(datahistlist[0],MCHistSum,1,1)
            elif len(DataFiles)==2: RatioHist.Divide(datahistlist[0],datahistlist[1],1,1)
            elif len(MCFiles)==2: RatioHist.Divide(mchistlist[0],MCHistSum,1,1)
            RatioHist.SetMinimum(0)
            RatioHist.SetMaximum(2)
            RatioHist.Draw("")
            MyLine = TLine(RatioHist.GetXaxis().GetXmin(),1,RatioHist.GetXaxis().GetXmax(),1)
            MyLine.Draw("")
            can.SaveAs("ZMuMuMVAStack_"+hist+label+"_Ratio.png")
        stack.Clear()
        leg.Clear()
        gDirectory.Delete("*")
