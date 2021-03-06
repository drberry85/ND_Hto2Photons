print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
#gROOT.Macro("$HOME/rootlogon.C")
gROOT.ProcessLine('.L /afs/cern.ch/work/n/nancy/private/Higgs_Moriond/CMSSW_5_3_7_patch4/src/ND_Hto2Photons/TreeReaders/tools/style-Egamma.C')
setEgammaStyle()
#egammaGrid(True)
gStyle.SetOptStat(000000)

print "Setting Initial Parameters."
#can = TCanvas("Plots","Plots",850,650)
can = TCanvas("Plots","Plots")
#can.SetLogy()
colors = [3, 2, 5, 7, 6]
datacolors = [1,2,3,4,5,6,7]
#leg = TLegend(0.15, 0.7, 0.3, 0.8)
#leg.SetFillColor(0)
#leg.SetBorderSize(0)
#leg.SetTextSize(0.03)
myText=TLatex()
myText.SetNDC()
myText.SetTextSize(0.035)


HistogramNames=["pho_pt","pho_idmva_rescaled","pho_tmva_photonid_Ephot","pho_tmva_photonid_Eregr","pho_tmva_photonid_sigmaOverE","pho_tmva_photonid_etawidth","pho_tmva_photonid_phiwidth","pho_tmva_photonid_r9","pho_tmva_photonid_r9_zoom","pho_mitmvaR9Bin1","pho_mitmvaR9Bin2","pho_mitmvaR9Bin3", "pho_tmva_photonid_sceta","pho_tmva_photonid_sieie","pho_tmva_photonid_sieip","pho_tmva_photonid_s4ratio","pho_tmva_photonid_lambdaratio","pho_tmva_photonid_eventrho","pho_tmva_photonid_ESEffSigmaRR","pho_tmva_photonid_pfchargedisogood03","pho_tmva_photonid_pfchargedisobad03","pho_tmva_photonid_pfphotoniso03","pho_tmva_photonid_pfneutraliso03"]
Categories=["Barrel","Endcap"]

HistogramNames2=["pho_mitmvaEtaBin1","pho_mitmvaEtaBin2","pho_mitmvaEtaBin3","pho_mitmvaEtaBin4","pho_tmva_photonid_sieipEndcapPlus","pho_tmva_photonid_sieipEndcapMinus"]
HistogramNames3=["pho_idmva_rescaled"]
pwd = "../../"
DataFileNames = ["ZMuMu_Run2012ABCD_PhoPFPresel_HighPt.root"]
#DataFileNames = ["ZMuMu_Run2012ABCD_HighPt.root"]
#DataFileNames = ["ZMuMu_Run2012ABCD_NoGSFVeto.root"]
DataLegends = ["8TeV Data "]
#DataLegends = ["8TeV Data, RunABC"]
#DataFileNames = ["ZMuMu_Run2012ABC_PhoPFPresel_HighPt.root","ZMuMu_Run2012D_PhoPFPresel_HighPt_noEcalIso.root"]
#DataFileNames = ["ZMuMu_Run2012ABC_PhoPFPresel_HighPt.root","ZMuMu_Run2012D_PhoPFPresel_HighPt.root"]
#DataLegends = ["Run2012ABC, 11.8fb","Run2012D, 7.3fb scaled to 11.8fb"]
#DataLegends = ["Run2012ABC, 11.8fb","Run2012D, 7.3fb scaled to ABC"]
FileNames = []
FileNames = ["ZMuMu_ZToMuMu_PhoPFPresel_Corr2012_HighPt.root"]
#FileNames = ["ZMuMu_ZToMuMu_Corr2012_HighPt.root"]
#FileNames = ["ZMuMu_ZToMuMu_NoGSFVeto_Corr2012.root"]
#FileNames = ["ZMuMu_ZToMuMu_PhoPFPresel_Corr2012_idMVAScalingABCD_HighPt.root"]
Legends = ["Z#rightarrow#mu#mu#gamma MC"]
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
        can.SetCanvasSize(600,600)
        if MakeLogScale: can.SetLogy()
        if histbase=="pho_tmva_photonid_sceta" and cat=="Endcap": hist=histbase
        else: hist=histbase+cat
        stack = THStack(hist,"")
        mchistlist=[]
        datahistlist=[]
        MCHistSum = TH1F()
        Normalization=0
        for i in range(len(MCFiles)): Normalization+=MCFiles[i].Get(hist).Integral()
        print hist
        for i in range(len(MCFiles)):
            MCFiles[i].cd()
            mchistlist.append(MCFiles[i].Get(hist))
            mchistlist[i].Scale(DataFiles[0].Get(hist).Integral()/Normalization)
            if i==0: MCHistSum=mchistlist[i]
            else: MCHistSum.Add(mchistlist[i])
            mchistlist[i].SetTitle("")
            if i==0: stack.SetTitle(hist+";"+mchistlist[0].GetXaxis().GetTitle()+";")
            mchistlist[i].SetLineWidth(0)
            mchistlist[i].SetLineColor(1)
            mchistlist[i].SetFillColor(colors[i])
            if hist.find("pfphotoniso") != -1: mchistlist[i].GetXaxis().SetTitle("Neutral EM isolation (GeV)")
            if hist.find("pfchargedisogood") != -1:
                mchistlist[i].GetXaxis().SetTitle("Charged hadron isolation (GeV) ")
                mchistlist[i].GetXaxis().SetTitleSize(0.01);
            if hist.find("pfchargedisobad") != -1:
                mchistlist[i].GetXaxis().SetTitleSize(0.002);
                mchistlist[i].GetXaxis().SetTitle("Max charged hadron isolation (GeV)")
            if hist.find("sieip") != -1: mchistlist[i].GetXaxis().SetTitle("Cov_{i#etai#phi}")
            if hist.find("etawidth") != -1: mchistlist[i].GetXaxis().SetTitle("#eta_{width}")
            if hist.find("phiwidth") != -1: mchistlist[i].GetXaxis().SetTitle("#phi_{width}")
            if hist.find("pho_pt") != -1:
                mchistlist[i].GetXaxis().SetTitle("P_{T} (GeV)")
            if hist.find("pho_idmva_rescaled") != -1: mchistlist[i].GetXaxis().SetTitle("Photon ID MVA")
            if hist.find("pho_tmva_photonid_sieipEndcap") != -1:
                leg = TLegend(0.15, 0.78, 0.3, 0.88)
            elif hist.find("sieieEndcap") != -1:
                leg = TLegend(0.15, 0.78, 0.3, 0.88)
            elif hist.find("etawidth") != -1:
                leg = TLegend(0.58, 0.78, 0.8, 0.88)
            elif hist.find("phiwidth") != -1:
                leg = TLegend(0.58, 0.78, 0.8, 0.88)
            elif hist.find("photoniso") != -1:
                leg = TLegend(0.58, 0.78, 0.8, 0.88)
            elif hist.find("chargedisobad") != -1:
                leg = TLegend(0.58, 0.78, 0.8, 0.88)
            elif hist.find("rho") != -1:
                leg = TLegend(0.58, 0.78, 0.8, 0.88)
            elif hist.find("sigmaOverE") != -1:
                leg = TLegend(0.58, 0.78, 0.8, 0.88)
            elif hist.find("pho_pt") != -1:
                leg = TLegend(0.58, 0.78, 0.8, 0.88)
            else:
                leg = TLegend(0.15, 0.7, 0.3, 0.8)
                
                
            leg.SetFillColor(0)
            leg.SetBorderSize(0)
            leg.SetTextSize(0.04)
            stack.Add(mchistlist[i])
            leg.AddEntry(mchistlist[i],Legends[i])
            if len(MCFiles)>0:
                # bw= str(mchistlist[0].GetBinWidth(3))
                bw= "%4.3g" %mchistlist[0].GetBinWidth(3)
                if hist.find("etawidth") != -1:  bw= "%6.4g" %mchistlist[0].GetBinWidth(3)
                if hist.find("sieie") != -1:  bw= "%5.4g" %mchistlist[0].GetBinWidth(3)
                if hist.find("sieip") != -1:  bw= "%7.6g" %mchistlist[0].GetBinWidth(3)
                stack.SetTitle(";"+mchistlist[0].GetXaxis().GetTitle()+";Events/"+bw)
                if hist.find("pho_pt") != -1:
                    bw= "%1.0g" %mchistlist[0].GetBinWidth(3)
                    stack.SetTitle(";"+mchistlist[0].GetXaxis().GetTitle()+";Events/"+bw+" GeV")
                if hist.find("pfchargedisogood") != -1:
                    bw= "%2.1g" %mchistlist[0].GetBinWidth(3)
                    stack.SetTitle(";"+mchistlist[0].GetXaxis().GetTitle()+";Events/"+bw+" GeV")
                if hist.find("pfchargedisobad") != -1:
                    bw= "%3.2g" %mchistlist[0].GetBinWidth(3)
                    stack.SetTitle(";"+mchistlist[0].GetXaxis().GetTitle()+";Events/"+bw+" GeV")
                if hist.find("photoniso") != -1:
                    bw= "%2.1g" %mchistlist[0].GetBinWidth(3)
                    stack.SetTitle(";"+mchistlist[0].GetXaxis().GetTitle()+";Events/"+bw+" GeV")
                stack.SetMaximum(1.3*stack.GetMaximum())
                stack.Draw('hist')
                if hist.find("pho_pt") != -1:
                    stack.GetXaxis().SetRangeUser(20,100)
                if hist.find("etawidth") != -1:
                    stack.GetXaxis().SetNdivisions(508);
                Normalization=0
                for i in range(len(DataFiles)):
                    if i==1: Normalization+=DataFiles[i].Get(hist).Integral()
                    for i in range(len(DataFiles)):
                        datahistlist.append(DataFiles[i].Get(hist))
                        #            datahistlist[i].Scale(DataFiles[0].Get(hist).Integral()/Normalization)
                        datahistlist[i].SetMarkerStyle(20)
                        datahistlist[i].SetMarkerColor(datacolors[i])
                        leg.AddEntry(datahistlist[i],DataLegends[i])
                        #if DataLegends[i].find("Run2012D")!=-1: datahistlist[i].Scale(11.8/7.3)
                        if i==0: datahistlist[0].SetMaximum(datahistlist[0].GetMaximum()*1.3)
                        if i==0 and len(MCFiles)==0: datahistlist[i].Draw("e")
                        else: datahistlist[i].Draw("esame")
                        leg.Draw()
                        myText.DrawLatex(0.15,0.92,"CMS preliminary, #sqrt{s} = 8 TeV,  L = 19.6 fb^{-1} ");
                        myText.SetTextSize(0.04)
                        if hist.find("rescaledBarrel") != -1:
                            myText.DrawLatex(0.17,0.65,"Barrel");
                        elif hist.find("rescaledEndcap") != -1:
                            myText.DrawLatex(0.17,0.65,"Endcap");
                        elif hist.find("good03Barrel") != -1:
                            myText.DrawLatex(0.17,0.65,"Barrel");
                        elif hist.find("good03Endcap") != -1:
                            myText.DrawLatex(0.17,0.65,"Endcap");
                        elif hist.find("good03Barrel") != -1:
                            myText.DrawLatex(0.17,0.65,"Barrel");
                        elif hist.find("good03Endcap") != -1:
                            myText.DrawLatex(0.17,0.65,"Endcap");
                        elif hist.find("r9_zoomBarrel") != -1:
                            myText.DrawLatex(0.17,0.65,"Barrel");
                        elif hist.find("r9_zoomEndcap") != -1:
                            myText.DrawLatex(0.17,0.65,"Endcap");
                        elif hist.find("r9Barrel") != -1:
                            myText.DrawLatex(0.17,0.65,"Barrel");
                        elif hist.find("r9Endcap") != -1:
                            myText.DrawLatex(0.17,0.65,"Endcap");
                        elif hist.find("s4ratioBarrel") != -1:
                            myText.DrawLatex(0.17,0.65,"Barrel");
                        elif hist.find("s4ratioEndcap") != -1:
                            myText.DrawLatex(0.17,0.65,"Endcap");
                        elif hist.find("sieieBarrel") != -1:
                            myText.DrawLatex(0.17,0.65,"Barrel");
                        elif hist.find("sieieEndcap") != -1:
                            myText.DrawLatex(0.17,0.65,"Endcap");
                        elif hist.find("sieipBarrel") != -1:
                            myText.DrawLatex(0.17,0.65,"Barrel");
                        elif hist.find("sieipEndcap") != -1:
                            myText.DrawLatex(0.17,0.65,"Endcap");
                        elif hist.find("pho_tmva_photonid_sieipEndcap") != -1:
                            myText.DrawLatex(0.17,0.65,"Endcap");
                        elif hist.find("sieieEndcap") != -1:
                            myText.DrawLatex(0.17,0.65,"Endcap");
                        elif hist.find("etawidthBarrel") != -1:
                            myText.DrawLatex(0.60,0.72,"Barrel");
                        elif hist.find("etawidthEndcap") != -1:
                            myText.DrawLatex(0.60,0.72,"Endcap");
                        elif hist.find("phiwidthBarrel") != -1:
                            myText.DrawLatex(0.60,0.72,"Barrel");
                        elif hist.find("phiwidthEndcap") != -1:
                            myText.DrawLatex(0.60,0.72,"Endcap");
                        elif hist.find("photonisoBarrel") != -1:
                            myText.DrawLatex(0.60,0.72,"Barrel");
                        elif hist.find("photonisoEndcap") != -1:
                            myText.DrawLatex(0.60,0.72,"Endcap");
                        elif hist.find("chargedisobad03Barrel") != -1:
                            myText.DrawLatex(0.60,0.72,"Barrel");
                        elif hist.find("chargedisobad03Endcap") != -1:
                            myText.DrawLatex(0.60,0.72,"Endcap");
                        elif hist.find("rhoBarrel") != -1:
                            myText.DrawLatex(0.60,0.72,"Barrel");
                        elif hist.find("rhoEndcap") != -1:
                            myText.DrawLatex(0.60,0.72,"Endcap");
                        elif hist.find("sigmaOverEBarrel") != -1:
                            myText.DrawLatex(0.60,0.72,"Barrel");
                        elif hist.find("sigmaOverEEndcap") != -1:
                            myText.DrawLatex(0.60,0.72,"Endcap");
                        elif hist.find("pho_ptBarrel") != -1:
                            myText.DrawLatex(0.60,0.72,"Barrel");
                        elif hist.find("pho_ptEndcap") != -1:
                            myText.DrawLatex(0.60,0.72,"Endcap");
                        can.SaveAs("ZMuMuMVAStack_"+hist+label+".png")
                        can.SaveAs("ZMuMuMVAStack_"+hist+label+".pdf")
                        can.SaveAs(hist+label+".C")
                        if len(DataFiles)+len(MCFiles)==2:
                            can.SetCanvasSize(600,300)
                            can.SetLogy(0)
                            if len(DataFiles)==0: RatioHist = mchistlist[0].Clone("")
                            else: RatioHist = datahistlist[0].Clone("")
                            if len(DataFiles)==1: RatioHist.Divide(datahistlist[0],MCHistSum,1,1)
                            elif len(DataFiles)==2: RatioHist.Divide(datahistlist[0],datahistlist[1],1,1)
                            elif len(MCFiles)==2: RatioHist.Divide(mchistlist[0],MCHistSum,1,1)
                            if hist.find("pfphotoniso") != -1: RatioHist.GetXaxis().SetTitle("Photon isolation")
                            if hist.find("pfchargedisogood") != -1: RatioHist.GetXaxis().SetTitle("Hadron isolation (right vertex)")
                            if hist.find("pfchargedisobad") != -1: RatioHist.GetXaxis().SetTitle("Hadron isolation (wrong vertex)")
                            RatioHist.SetMinimum(0)
                            RatioHist.SetMaximum(2)
                            RatioHist.SetTitle("")
                            RatioHist.GetYaxis().SetTitle("Ratio");
                            RatioHist.Draw("")
                            MyLine = TLine(RatioHist.GetXaxis().GetXmin(),1,RatioHist.GetXaxis().GetXmax(),1)
                            MyLine.Draw("")
                            can.SaveAs("ZMuMuMVAStack_"+hist+label+"_Ratio.png")
                            can.SaveAs("ZMuMuMVAStack_"+hist+label+"_Ratio.pdf")
                            can.SaveAs(hist+label+"_Ratio.C")
                            stack.Clear()
                            leg.Clear()
                            


#for histbase in HistogramNames2:
#        can.SetCanvasSize(600,600)
#        if MakeLogScale: can.SetLogy()
#        hist=histbase
#        stack = THStack(hist,"")
#        mchistlist=[]
#        datahistlist=[]
#        MCHistSum = TH1F()
#        Normalization=0
#        for i in range(len(MCFiles)): Normalization+=MCFiles[i].Get(hist).Integral()
#        for i in range(len(MCFiles)):
#            MCFiles[i].cd()
#            mchistlist.append(MCFiles[i].Get(hist))
#            mchistlist[i].Scale(DataFiles[0].Get(hist).Integral()/Normalization)
#            if i==0: MCHistSum=mchistlist[i]
#            else: MCHistSum.Add(mchistlist[i])
#            mchistlist[i].SetTitle("")
#            if i==0: stack.SetTitle(hist+";"+mchistlist[0].GetXaxis().GetTitle()+";")
#            mchistlist[i].SetLineWidth(0)
#            mchistlist[i].SetLineColor(colors[i])
#            mchistlist[i].SetFillColor(colors[i])
#            stack.Add(mchistlist[i])
#            leg.AddEntry(mchistlist[i],Legends[i])
#        if len(MCFiles)>0:
#            stack.SetTitle(";"+mchistlist[0].GetXaxis().GetTitle()+";Events in Data")
#            stack.SetMaximum(1.3*stack.GetMaximum())
#            stack.Draw('hist')
#            Normalization=0
#        for i in range(len(DataFiles)):
#            if i==1: Normalization+=DataFiles[i].Get(hist).Integral()
#        for i in range(len(DataFiles)):
#            datahistlist.append(DataFiles[i].Get(hist))
##            datahistlist[i].Scale(DataFiles[0].Get(hist).Integral()/Normalization)
#            datahistlist[i].SetMarkerStyle(20)
#            datahistlist[i].SetMarkerColor(datacolors[i])
#            leg.AddEntry(datahistlist[i],DataLegends[i])
#            #if DataLegends[i].find("Run2012D")!=-1: datahistlist[i].Scale(11.8/7.3)
#            if i==0: datahistlist[0].SetMaximum(datahistlist[0].GetMaximum()*1.3)
#            if i==0 and len(MCFiles)==0: datahistlist[i].Draw("e")
#            else: datahistlist[i].Draw("esame")
#            leg.Draw()
#        can.SaveAs("ZMuMuMVAStack_"+hist+label+".png")
#        if len(DataFiles)+len(MCFiles)==2:
#            can.SetCanvasSize(600,300)
#            can.SetLogy(0)
#            if len(DataFiles)==0: RatioHist = mchistlist[0].Clone("")
#            else: RatioHist = datahistlist[0].Clone("")
#            if len(DataFiles)==1: RatioHist.Divide(datahistlist[0],MCHistSum,1,1)
#            elif len(DataFiles)==2: RatioHist.Divide(datahistlist[0],datahistlist[1],1,1)
#            elif len(MCFiles)==2: RatioHist.Divide(mchistlist[0],MCHistSum,1,1)
#            RatioHist.SetMinimum(0)
#            RatioHist.SetMaximum(2)
#            RatioHist.Draw("")
#            MyLine = TLine(RatioHist.GetXaxis().GetXmin(),1,RatioHist.GetXaxis().GetXmax(),1)
#            MyLine.Draw("")
#            can.SaveAs("ZMuMuMVAStack_"+hist+label+"_Ratio.png")
#        stack.Clear()
#        leg.Clear()



        
        gDirectory.Delete("*")
