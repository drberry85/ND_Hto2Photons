print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
import array
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

print "Setting Initial Parameters."
can = TCanvas("Plots","Plots",850,650)
can.SetGrid()
leg = TLegend(0.7, 0.8, 0.9, 0.9)
leg.SetFillColor(0)
leg.SetBorderSize(1)
#leg.SetNColumns(3)
#HistogramNames=["ConvdZEta","SuperdZEta","MixdZEta"]
HistogramNames=["MixdZEta"]
legendlist=[["Data RunA+B"],["#gamma+Jet MC"]]
#legendlist=[["1cm Data","5mm Data","3mm Data"],["1cm #gamma+Jet MC","5mm #gamma+Jet MC","3mm #gamma+Jet MC"],["1cm 120GeV Higgs","5mm 120GeV Higgs","3mm 120GeV Higgs"]]
#legendlist=[["1cm Signal 42X","5mm Signal 42X","3mm Signal 42X"],["1cm Signal 44X","5mm Signal 44X","3mm Signal 44X"]]
OutFileNames=["MixEtaEffData.png"]
#OutFileNames=["SuperClusterEtaEff44X.png","ConvEtaEff44X.png","MixEtaEff44X.png"]
filenamelist=["done/Vertex_Data.root","better/Vertex_PhotonPlusJetMC.root"]
#filenamelist=["Vertex_PhotonPlusJetData.root","Vertex_PhotonPlusJetMC.root","Vertex_120GeV_SimVertex.root"]
#filenamelist=["Vertex_120GeV_SimVertex.root","Vertex_130GeV_SimVertex.root"]
pwd = "/data/ndpc2/c/HiggsGammaGamma/PhotonPlusJet/CMSSW_4_2_3/src/ND_Hto2Photons/TreeReaders/"
lowbins=[41]
upperbins=[60]
#lowbins=[41,46,48]
#upperbins=[60,55,53]
graphcolor=[1, 2, 4, 7, 3, 6, 5, 9, 8]
MarkerStyles=[20,20]
LineStyles=[1,2]
#MarkerStyles=[20,21,22,23,24]
#LineStyles=[1,2,3,4]

for histname,outfile in zip(HistogramNames,OutFileNames):
	multigraph = TMultiGraph()
	for filename,legend,linestyle,markerstyle,color2 in zip(filenamelist,legendlist,LineStyles,MarkerStyles,graphcolor):
		file=TFile(pwd+filename)
		hist=file.Get(histname)
		for Min,Max,legendname,color in zip(lowbins,upperbins,legend,graphcolor):
			XValue=[]
			XError=[]
			HistEff=[]
			HistError=[]
			for i in range(1,hist.GetNbinsX()+1):
				XValue.append(float(i)/10-0.05)
				XError.append(0.05)
				if hist.Integral(i,i,Min,Max)!=0:
					HistEff.append(hist.Integral(i,i,Min,Max)/hist.Integral(i,i,1,100)*100)
					HistError.append(sqrt(hist.Integral(i,i,Min,Max))/hist.Integral(i,i,1,100)*100)
				else:
					HistEff.append(0)
					HistError.append(0)
				#print "XValue: %f XError: %f HistEff: %f HistError: %f" %(XValue[i-1],XError[i-1],HistEff[i-1],HistError[i-1])
			points = array.array('f',XValue)
			pointserror = array.array('f',XError)
			data = array.array('f',HistEff)
			dataerror = array.array('f',HistError)
			graph=TGraphErrors(hist.GetNbinsX(),points,data,pointserror,dataerror)
			graph.SetLineColor(color2)
			graph.SetMarkerStyle(markerstyle)
			graph.SetMarkerSize(1)
			graph.SetLineStyle(linestyle)
			graph.SetMarkerColor(color2)
			graph.SetFillColor(kWhite)
			leg.AddEntry(graph,legendname)
			multigraph.Add(graph)
	multigraph.SetNameTitle("multigraph",";#eta;Efficiency (%)");
	multigraph.SetMinimum(75)
	multigraph.SetMaximum(100)
	multigraph.Draw('AP')
	leg.Draw()
	can.SaveAs(outfile)
	can.Clear()
	leg.Clear()
	multigraph.Delete("*")

print "Done"
