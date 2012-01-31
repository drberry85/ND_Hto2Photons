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
leg = TLegend(0.6, 0.6, 0.9, 0.9)
leg.SetFillColor(0)
leg.SetBorderSize(1)
pwd = "/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/"
hists1 = ["PixelBarrelConvdZ","PixelFwdConvdZ","TECSuperdZ","TIBSuperdZ","TIDConvdZ","TOBSuperdZ","ConvdZ"]
#hists1 = ["PixelBarrelConvdZRes","PixelFwdConvdZRes","TECSuperdZRes","TIBSuperdZRes","TIDConvdZRes","TOBSuperdZRes"]
#hists1 = ["PixelBarrelConvdZ","PixelFwdConvdZ","TECConvdZ","TIBConvdZ","TIDConvdZ","TOBConvdZ","ConvdZ"]
#hists2 = ["PixelBarrelSuperdZ","PixelFwdSuperdZ","TECSuperdZ","TIBSuperdZ","TIDSuperdZ","TOBSuperdZ","SuperdZ"]
titles = ["Pixel Barrel","Pixel Forward","Tracker EndCap","Tracker Inner Barrel","Tracker Inner Disk","Tracker Outer Barrel","All Tracker Regions"]
Data = TFile(pwd+"Vertex_Data.root")
MC = TFile(pwd+"Vertex_PhotonPlusJetMC.root")

for hist1,title in zip(hists1,titles):
	DataHist1=Data.Get(hist1)
	MCHist1=MC.Get(hist1)
	#MCHist1.Sumw2()
	#DataHist2=Data.Get(hist2)
	#MCHist2=MC.Get(hist2)
	#MCHist2.Sumw2()
	DataHist1.Sumw2()
	DataHist1.Scale(1/DataHist1.Integral())
	DataHist1.SetMarkerColor(kBlue)
	#DataHist2.Sumw2()
	#DataHist2.Scale(1/DataHist2.Integral())
	#DataHist2.SetMarkerColor(kGreen)
	MCHist1.Scale(1/MCHist1.Integral())
	MCHist1.SetLineColor(kBlue)
	#MCHist1.SetLineWidth(2)
	MCHist1.SetFillStyle(3004)
	MCHist1.SetFillColor(kBlue)
	#MCHist2.Scale(1/MCHist2.Integral())
	#MCHist2.SetLineColor(kGreen)
	#MCHist2.SetLineWidth(2)
	#MCHist2.SetFillStyle(3005)
	#MCHist2.SetFillColor(kGreen)
	leg.AddEntry(MCHist1,"#gamma+Jet Monte Carlo")
	#leg.AddEntry(MCHist2,"#splitline{Regular Conversions }{+ SuperCluster}")
	#MCHist2.SetMaximum(max(MCHist1.GetMaximum(),MCHist2.GetMaximum(),DataHist1.GetMaximum(),DataHist2.GetMaximum())*1.1)
	#MCHist2.SetTitle("")
	#MCHist2.Draw("")
	MCHist1.Draw("")
	DataHist1.SetMarkerStyle(20)
	#DataHist2.SetMarkerStyle(20)
	#MCHist1.GetXaxis().SetRangeUser(-5,5)
	#DataHist1.GetXaxis().SetRangeUser(-5,5)
	leg.AddEntry(DataHist1,"Run 2011A+B")
	#leg.AddEntry(DataHist2,"#splitline{Single Leg Conversions }{+ SuperCluster}")
	DataHist1.Draw("esame")
	#DataHist2.Draw("esame")
	leg.Draw()
	can.SaveAs(hist1+".png")
	leg.Clear()

# 	MCHist1Check=MC.Get(hist1+"sim")
# 	MCHist1Check.Scale(1/MCHist1Check.Integral())
# 	MCHist1Check.SetLineColor(kRed)
# 	MCHist1Check.SetLineWidth(2)
# 	MCHist1.SetLineColor(kBlack)
# 	MCHist1.Draw("")
# 	MCHist1Check.Draw("same")
# 	leg.AddEntry(MCHist1,"PJet MC Conversion")
# 	leg.AddEntry(MCHist1Check,"PJet MC SimVertex")
# 	leg.Draw()
# 	can.SaveAs(hist1+"Check.png")
# 	leg.Clear()
	
# 	MCHist2Check=MC.Get(hist2+"sim")
# 	MCHist2Check.Scale(1/MCHist2Check.Integral())
# 	MCHist2Check.SetLineColor(kRed)
# 	MCHist2Check.SetLineWidth(2)
# 	MCHist2.SetLineColor(kBlack)
# 	MCHist2.Draw("")
# 	MCHist2Check.Draw("same")
# 	leg.AddEntry(MCHist2,"PJet MC SuperCluster")
# 	leg.AddEntry(MCHist2Check,"PJet MC SimVertex")
# 	leg.Draw()
# 	can.SaveAs(hist2+"Check.png")
# 	leg.Clear()
	
	gDirectory.Delete("*")

print "Done!"
