print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
import array
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

print "Setting Initial Parameters."
can = TCanvas("Plots","Plots",850,650)
can.SetGrid()
can.SetFillColor(0)
leg = TLegend(0.37, 0.13, 0.87, 0.28)
leg.SetFillColor(0)
leg.SetBorderSize(1)
leg.SetNColumns(3)
SuperClusterNames=["PixelBarrelSuperdZ","TIBSuperdZ","TOBSuperdZ","PixelFwdSuperdZ","TIDSuperdZ","TECSuperdZ"]
ConvNames=["PixelBarrelConvdZ","TIBConvdZ","TOBConvdZ","PixelFwdConvdZ","TIDConvdZ","TECConvdZ"]
MixedConversion=["PixelBarrelConvdZ","TIBSuperdZ","TOBSuperdZ","PixelFwdConvdZ","TIDConvdZ","TECSuperdZ"]
#SuperClusterNames.append("dZcheck")
#ConvNames.append("dZcheck")
#MixedConversion.append("dZcheck")
HistogramNames=[MixedConversion,MixedConversion]
#HistogramNames=[MixedConversion]
#legendlist=[["Pixel Barrel","TIB","TOB","Pixel Forward","TID","TEC"],["Pixel Barrel","TIB","TOB","Pixel Forward","TID","TEC"],["Pixel Barrel","TIB","TOB","Pixel Forward","TID","TEC"]]
legendlist=["Pixel Barrel","TIB","TOB","Pixel Forward","TID","TEC"]
#legendlist=[["Pix Barrel Data","TIB Data","TOB Data","Pix FWD Data","TID Data","TEC Data","40GeV Jet"],["Pix Barrel #gamma+Jet MC","TIB #gamma+Jet MC","TOB #gamma+Jet MC","Pix FWD #gamma+Jet MC","TID #gamma+Jet MC","TEC #gamma+Jet MC","40GeV Jet #gamma+Jet MC"]]
#legendlist=[["Pix Barrel Data","TIB Data","TOB Data","Pix FWD Data","TID Data","TEC Data","40GeV Jet"],["Pix Barrel SingleLeg","TIB SingleLeg","TOB SingleLeg","Pix FWD SingleLeg","TID SingleLeg","TEC SingleLeg","40GeV Jet SingleLeg"]]
OutFileNames=["MixedEffData.png","MixedEffPJet.png"]
#OutFileNames=["SuperEff_Higgs125GeV_52X.png","SuperEff_Higgs125GeV_52X_DoubleLeg.png","SuperEff_Higgs125GeV_52X_SingleLeg.png","ConvEff_Higgs125GeV_52X.png","ConvEff_Higgs125GeV_52X_DoubleLeg.png","ConvEff_Higgs125GeV_52X_SingleLeg.png"]
#OutFileNames=["MixedEffPJetSummer11.png"]
#filenamelist=["Vertex_PhotonPlusJetMC_Summer11.root"]
filenamelist=["Vertex_Data.root","Vertex_PJet.root"]
#filenamelist=["Vertex_Higgs_125GeV52X_SimVertex.root","Vertex_Higgs_125GeV52X_SimVertex_DoubleLeg.root","Vertex_Higgs_125GeV52X_SimVertex_SingleLeg.root","Vertex_Higgs_125GeV52X_SimVertex.root","Vertex_Higgs_125GeV52X_SimVertex_DoubleLeg.root","Vertex_Higgs_125GeV52X_SimVertex_SingleLeg.root"]
pwd = "/data/ndpc2/b/drberry/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/"
graphcolor=[862,870,433,616,880,619,1]
markers=[20,21,22,23,24,25,26]

for histlist,outfile,filename in zip(HistogramNames,OutFileNames,filenamelist):
	multigraph = TMultiGraph()
	file=TFile(pwd+filename)
	for histname,color,marstyle,legendname in zip(histlist,graphcolor,markers,legendlist):
		hist = file.Get(histname)
		HistIntegral = hist.Integral()+hist.GetBinContent(0)+hist.GetBinContent(101)
		print "Looking at Histogram %s" %histname
		print "%s Total: %i" %(histname,int(HistIntegral))
		print "%s 10mm: %i" %(histname,int(hist.Integral(41,60)))
		print "%s 5mm: %i" %(histname,int(hist.Integral(46,55)))
		print "%s 3mm: %i" %(histname,int(hist.Integral(48,53)))
		points = array.array('f',[10,5,3])
		pointserror = array.array('f',[0,0,0])
		data = array.array('f',[hist.Integral(41,60)/HistIntegral*100, hist.Integral(46,55)/HistIntegral*100, hist.Integral(48,53)/HistIntegral*100])
		dataerror = array.array('f',[sqrt(hist.Integral(41,60))/HistIntegral*100, sqrt(hist.Integral(46,55))/HistIntegral*100, sqrt(hist.Integral(48,53))/HistIntegral*100])
		graph=TGraphErrors(3,points,data,pointserror,dataerror)
		graph.SetLineColor(color)
		graph.SetMarkerStyle(marstyle)
		graph.SetMarkerSize(1.5)
		graph.SetLineWidth(2)
		#if (legendname.find("42X")==-1 and legendname.find("Data")==-1): graph.SetLineStyle(2)
		graph.SetMarkerColor(color)
		graph.SetFillColor(kWhite)
		leg.AddEntry(graph,legendname)
		multigraph.Add(graph)
	multigraph.SetNameTitle("multigraph",";#DeltaZ Cut (mm);Efficiency (%)");
	multigraph.SetMinimum(0)
	multigraph.SetMaximum(100)
	multigraph.Draw('ALP')
	leg.Draw()
	can.SaveAs(outfile)
	can.Clear()
	leg.Clear()
	multigraph.Delete("*")

print "Done"
