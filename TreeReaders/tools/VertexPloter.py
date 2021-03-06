print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)
gStyle.SetCanvasBorderMode(0);
gStyle.SetCanvasColor(kWhite);

print "Setting Initial Parameters."
can = TCanvas("Plots","Plots",850,650)
#can.SetLogy()
file=[]
leg = TLegend(0.6, 0.75, 0.9, 0.9)
leg.SetFillColor(0)
leg.SetBorderSize(1)
pwd = "../"
#hists1 = ["PixelBarrelConvdZEff","PixelFwdConvdZEff","TECConvdZEff","TIBConvdZEff","TIDConvdZEff","TOBConvdZEff","ConvdZ"]
hists1 = ["PixelBarrelConvdZEff","PixelFwdConvdZEff","TECSuperdZEff","TIBSuperdZEff","TIDConvdZEff","TOBSuperdZEff"]
#hists2 = ["PixelBarrelSuperdZEff","PixelFwdSuperdZEff","TECSuperdZEff","TIBSuperdZEff","TIDSuperdZEff","TOBSuperdZEff","SuperdZ"]
#hists1 = ["PixelBarrelRishiTrackdZEff","PixelFwdRishiTrackdZEff","TECRishiTrackdZEff","TIBRishiTrackdZEff","TIDRishiTrackdZEff","TOBRishiTrackdZEff","RishiTrackdZ"]
#hists2 = ["PixelBarrelRishiSuperdZEff","PixelFwdRishiSuperdZEff","TECRishiSuperdZEff","TIBRishiSuperdZEff","TIDRishiSuperdZEff","TOBRishiSuperdZEff","RishiSuperdZ"]
#hists2 = ["PixelBarrelRishiCombineddZEff","PixelFwdRishiCombineddZEff","TECRishiCombineddZEff","TIBRishiCombineddZEff","TIDRishiCombineddZEff","TOBRishiCombineddZEff","RishiCombineddZ"]
Data = TFile(pwd+"Vertex_Data.root")
MC = TFile(pwd+"Vertex_Background.root")

for hist1 in hists1:
    DataHist1=Data.Get(hist1)
    MCHist1=MC.Get(hist1)
    #DataHist2=Data.Get(hist2)
    #MCHist2=MC.Get(hist2)
    DataHist1.Scale(1/DataHist1.Integral())
    DataHist1.SetMarkerColor(kBlue)
    #DataHist2.Scale(1/DataHist2.Integral())
    #DataHist2.SetMarkerColor(kGreen)
    MCHist1.Scale(1/MCHist1.Integral())
    MCHist1.SetLineColor(kBlue)
    MCHist1.SetLineWidth(2)
    MCHist1.SetFillStyle(3004)
    MCHist1.SetFillColor(kBlue)
    #MCHist2.Scale(1/MCHist2.Integral())
    #MCHist2.SetLineColor(kGreen)
    #MCHist2.SetLineWidth(2)
    #MCHist2.SetFillStyle(3005)
    #MCHist2.SetFillColor(kGreen)
    leg.AddEntry(MCHist1,"#gamma+Jet and Background Monte Carlo")
    #leg.AddEntry(MCHist2,"Supercluster Higgs MC")
    #MCHist2.SetMaximum(max(DataHist1.GetMaximum(),DataHist2.GetMaximum(),MCHist1.GetMaximum(),MCHist2.GetMaximum())*1.2)
    #MCHist2.SetMaximum(max(MCHist1.GetMaximum(),MCHist2.GetMaximum())*1.2)
    #MCHist2.SetTitle("")
    #MCHist2.Draw("hist")
    #MCHist1.Draw("hist same")
    MCHist1.SetTitle("")
    MCHist1.Draw("hist")
    DataHist1.SetMarkerStyle(20)
    #DataHist2.SetMarkerStyle(20)
    DataHist1.Draw("esame")
    MCHist1.SetMaximum(1.1*max(MCHist1.GetMaximum(),DataHist1.GetMaximum()))
    leg.AddEntry(DataHist1,"Data Run2012A+B+C+D")
    #DataHist2.Draw("esame")
    #leg.AddEntry(DataHist2,"Supercluster Run2012 Data")
    leg.Draw()
    can.SaveAs("Vertex_"+hist1+"_Note.png")
    leg.Clear()

#   MCHist1Check=MC.Get(hist1+"sim")
#   MCHist1Check.Scale(1/MCHist1Check.Integral())
#   MCHist1Check.SetLineColor(kRed)
#   MCHist1Check.SetLineWidth(2)
#   MCHist1.SetLineColor(kBlack)
#   MCHist1.Draw("")
#   MCHist1Check.Draw("same")
#   leg.AddEntry(MCHist1,"PJet MC Conversion")
#   leg.AddEntry(MCHist1Check,"PJet MC SimVertex")
#   leg.Draw()
#   can.SaveAs(hist1+"Check.png")
#   leg.Clear()
    
#   MCHist2Check=MC.Get(hist2+"sim")
#   MCHist2Check.Scale(1/MCHist2Check.Integral())
#   MCHist2Check.SetLineColor(kRed)
#   MCHist2Check.SetLineWidth(2)
#   MCHist2.SetLineColor(kBlack)
#   MCHist2.Draw("")
#   MCHist2Check.Draw("same")
#   leg.AddEntry(MCHist2,"PJet MC SuperCluster")
#   leg.AddEntry(MCHist2Check,"PJet MC SimVertex")
#   leg.Draw()
#   can.SaveAs(hist2+"Check.png")
#   leg.Clear()
    
    gDirectory.Delete("*")

print "Done!"
