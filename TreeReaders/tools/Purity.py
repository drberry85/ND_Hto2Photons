print "Loading Root..."

from ROOT import *
import pdb
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

print "Setting Initial Parameters."
pwd = "/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_8_5_patch3/src/ND_Hto2Photons/TreeReaders/"
SignalFiles = ["HiggsAnalysis110GeV.root","HiggsAnalysis115GeV.root","HiggsAnalysis120GeV.root","HiggsAnalysis130GeV.root","HiggsAnalysis140GeV.root","HiggsAnalysis150GeV.root"]
BackgroundFiles = ["PhotonPlusJet.root","DoubleEMEnriched.root","QCDBCtoE.root","Box.root"]
hists = [["h_mass_2gammaAllEB","h_mass_2gammaAllEE","h_mass_2gammaSelEB","h_mass_2gammaSelEE","h_mass_2gammaMatchedEB","h_mass_2gammaMatchedEE"]]
hists.append(["h_mass_2gammaGoldenAllEB","h_mass_2gammaGoldenAllEE","h_mass_2gammaGoldenSelEB","h_mass_2gammaGoldenSelEE","h_mass_2gammaGoldenMatchedEB","h_mass_2gammaGoldenMatchedEE"])
hists.append(["h_mass_2gamma1goodconvAllEB","h_mass_2gamma1goodconvAllEE","h_mass_2gamma1goodconvSelEB","h_mass_2gamma1goodconvSelEE","h_mass_2gamma1goodconvMatchedEB","h_mass_2gamma1goodconvMatchedEE"])
hists.append(["h_mass_2gamma1poorconvAllEB","h_mass_2gamma1poorconvAllEE","h_mass_2gamma1poorconvSelEB","h_mass_2gamma1poorconvSelEE","h_mass_2gamma1poorconvMatchedEB","h_mass_2gamma1poorconvMatchedEE"])
hists.append(["h_mass_2gamma2convAllEB","h_mass_2gamma2convAllEE","h_mass_2gamma2convSelEB","h_mass_2gamma2convSelEE","h_mass_2gamma2convMatchedEB","h_mass_2gamma2convMatchedEE"])
hists.append(["h_mass_2gammaleftoverAllEB","h_mass_2gammaleftoverAllEE","h_mass_2gammaleftoverSelEB","h_mass_2gammaleftoverSelEE","h_mass_2gammaleftoverMatchedEB","h_mass_2gammaleftoverMatchedEE"])
background = []
data = []

print "Loading Background"
for i in range (len(BackgroundFiles)):
	background.append(TFile(pwd+BackgroundFiles[i]))

for i in range(len(SignalFiles)):
	backgroundmatrix = []
	backgroundsum = []
	eventmatrix = []
	print "Opening file: %s" %(SignalFiles[i])
	data.append(TFile(pwd+SignalFiles[i]))
	HiggsMass = SignalFiles[i][SignalFiles[i].find("Analysis")+8:SignalFiles[i].find("GeV")]
	LowerBin=int(floor(float(HiggsMass)*0.98)-80)
	UpperBin=int(floor(float(HiggsMass)*1.02)+1-80)
	print "The Mass Range is: %d to %d" %(float(HiggsMass)*0.98,float(HiggsMass)*1.02)
	print "Mass Range is bin number %d to %d" %(LowerBin,UpperBin)
	for j in range(len(BackgroundFiles)):
		print "Looking at Background File: %s" %BackgroundFiles[j]
		background[j].cd()
		backgroundlist=[]
		for k in range(len(hists)):
			events=[]
			hist=[]
			for l in range(len(hists[k])):
				hist.append(gDirectory.Get(hists[k][l]))
				#import pdb; pdb.set_trace()
				events.append(hist[l].Integral(LowerBin,UpperBin))
				print "The number of events in %s is: %0.02f" %(hists[k][l],events[l])
			backgroundlist.append(events)
		backgroundmatrix.append(backgroundlist)
	data[i].cd()
	for j in range(len(hists)):
		backgroundsum.append([0.0]*len(hists[j]))
	for j in range(len(backgroundsum)):
		for k in range(len(backgroundsum[k])):
			for l in range(len(backgroundmatrix)):
				backgroundsum[j][k] = backgroundsum[j][k] + backgroundmatrix[l][j][k]
	for j in range(len(hists)):
		print "Looking at file %s: %s" %(SignalFiles[i],hists[j][0][hists[j][0].find("2"):hists[j][0].find("All")])
		events=[]
		hist=[]
		for k in range(len(hists[j])):
			hist.append(gDirectory.Get(hists[j][k]))
			events.append(hist[k].Integral(LowerBin,UpperBin))
			print "The number of events in %s is: %0.02f with %.02f background events so the purity is %.02f%%" %(hists[j][k],events[k],backgroundsum[j][k],100*events[k]/(events[k]+backgroundsum[j][k]))
		eventmatrix.append(events)
	print "\n"
	print "Results Table for %s" %SignalFiles[i]
	print "Selection\t\t%s\t%s\t\t%s\t\t%s\t\t\t%s" %(SignalFiles[i],BackgroundFiles[0],BackgroundFiles[1],BackgroundFiles[2],BackgroundFiles[3])
	for j in range(len(eventmatrix)):
		if (len(hists[j][0][8:hists[j][0].find("All")])<6):
			print "%s:\t\t\t%.03e (%.03e)\t\t%.03e (%.03e)\t\t%.03e (%.03e)\t\t%.03e (%.03e)\t\t%.03e (%.03e)" %(hists[j][0][8:hists[j][0].find("All")],eventmatrix[j][2],eventmatrix[j][3],backgroundmatrix[0][j][2],backgroundmatrix[0][j][3],backgroundmatrix[1][j][2],backgroundmatrix[1][j][3],backgroundmatrix[2][j][2],backgroundmatrix[2][j][3],backgroundmatrix[3][j][2],backgroundmatrix[3][j][3])
		else:
			print "%s:\t\t%.03e (%.03e)\t\t%.03e (%.03e)\t\t%.03e (%.03e)\t\t%.03e (%.03e)\t\t%.03e (%.03e)" %(hists[j][0][8:hists[j][0].find("All")],eventmatrix[j][2],eventmatrix[j][3],backgroundmatrix[0][j][2],backgroundmatrix[0][j][3],backgroundmatrix[1][j][2],backgroundmatrix[1][j][3],backgroundmatrix[2][j][2],backgroundmatrix[2][j][3],backgroundmatrix[3][j][2],backgroundmatrix[3][j][3])
	print "\n"
	print "Purity Results for %s" %SignalFiles[i]
	print "Purity Table:\t\tSignal Events\t\t\tBackground Events\t\tPurity"
	for j in range(len(eventmatrix)):
		if (len(hists[j][0][8:hists[j][0].find("All")])<6):
			print "%s:\t\t\t%.03e (%.03e)\t\t%.03e (%.03e)\t\t%.02f%% (%.02f%%)" %(hists[j][0][8:hists[j][0].find("All")],eventmatrix[j][2],eventmatrix[j][3],backgroundsum[j][2],backgroundsum[j][3],100*eventmatrix[j][2]/(backgroundsum[j][2]+eventmatrix[j][2]),100*eventmatrix[j][3]/(backgroundsum[j][3]+eventmatrix[j][3]))
		else:
			print "%s:\t\t%.03e (%.03e)\t\t%.03e (%.03e)\t\t%.02f%% (%.02f%%)" %(hists[j][0][8:hists[j][0].find("All")],eventmatrix[j][2],eventmatrix[j][3],backgroundsum[j][2],backgroundsum[j][3],100*eventmatrix[j][2]/(backgroundsum[j][2]+eventmatrix[j][2]),100*eventmatrix[j][3]/(backgroundsum[j][3]+eventmatrix[j][3]))

print "Done!"
#import pdb; pdb.set_trace()
