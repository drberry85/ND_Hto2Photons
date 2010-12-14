print "Loading Root..."

from ROOT import *
import pdb
from math import log
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

GlobalWeight = 1000.0
EventReference = [["HiggsAnalysis110GeV.root",GlobalWeight*0.001939*(20.493+1.4405+1.4421)],["HiggsAnalysis115GeV.root",GlobalWeight*0.002101*(18.735+1.3712+1.2524)],["HiggsAnalysis120GeV.root",GlobalWeight*0.002219*(17.173+1.3062+1.0921)],["HiggsAnalysis130GeV.root",GlobalWeight*0.002240*(14.579+1.1866+0.8395)],["HiggsAnalysis140GeV.root",GlobalWeight*0.001929*(12.525+1.0811+0.6539)],["HiggsAnalysis150GeV.root",GlobalWeight*0.001363*(10.863+0.9868+0.5155)]]

print "Setting Initial Parameters."
pwd = "/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_8_5_patch3/src/ND_Hto2Photons/TreeReaders/"
SignalFiles = ["HiggsAnalysis110GeV.root","HiggsAnalysis115GeV.root","HiggsAnalysis120GeV.root","HiggsAnalysis130GeV.root","HiggsAnalysis140GeV.root","HiggsAnalysis150GeV.root"]
#SignalFiles = ["HiggsAnalysis110GeV.root"]
BackgroundFiles = ["PhotonPlusJet.root","ReweightedDoubleEMEnriched.root","QCDBCtoE.root","Box.root"]
hists = [["h_mass_2gammaAllEB","h_mass_2gammaAllEE","h_mass_2gammaSelEB","h_mass_2gammaSelEE","h_mass_2gammaMatchedEB","h_mass_2gammaMatchedEE"]]
hists.append(["h_mass_2gammaGoldenAllEB","h_mass_2gammaGoldenAllEE","h_mass_2gammaGoldenSelEB","h_mass_2gammaGoldenSelEE","h_mass_2gammaGoldenMatchedEB","h_mass_2gammaGoldenMatchedEE"])
hists.append(["h_mass_2gamma1goodconvAllEB","h_mass_2gamma1goodconvAllEE","h_mass_2gamma1goodconvSelEB","h_mass_2gamma1goodconvSelEE","h_mass_2gamma1goodconvMatchedEB","h_mass_2gamma1goodconvMatchedEE"])
hists.append(["h_mass_2gamma1poorconvAllEB","h_mass_2gamma1poorconvAllEE","h_mass_2gamma1poorconvSelEB","h_mass_2gamma1poorconvSelEE","h_mass_2gamma1poorconvMatchedEB","h_mass_2gamma1poorconvMatchedEE"])
hists.append(["h_mass_2gamma2convAllEB","h_mass_2gamma2convAllEE","h_mass_2gamma2convSelEB","h_mass_2gamma2convSelEE","h_mass_2gamma2convMatchedEB","h_mass_2gamma2convMatchedEE"])
hists.append(["h_mass_2gammaleftoverAllEB","h_mass_2gammaleftoverAllEE","h_mass_2gammaleftoverSelEB","h_mass_2gammaleftoverSelEE","h_mass_2gammaleftoverMatchedEB","h_mass_2gammaleftoverMatchedEE"])
background = []
data = []

fit=TF1("fit","[0]*exp([1]*(x-100))",100,160)
fit.SetParNames("constant","coefficient");
fit.SetParLimits(1,-1,0)
fit.SetLineColor(4)

print "Loading Background"
for i in range (len(BackgroundFiles)):
	background.append(TFile(pwd+BackgroundFiles[i]))

for i in range(len(SignalFiles)):
	backgroundmatrix = []
	backgroundsum = []
	backgroundintegral = []
	backgroundintegralerror = []
	backgroundhist = []
	eventmatrix = []
	for j in range(len(hists)):
		backgroundhisttemp = []
		backgroundsum.append([0.0]*len(hists[j]))
		backgroundintegral.append([0.0]*len(hists[j]))
		backgroundintegralerror.append([0.0]*len(hists[j]))
		for k in range(len(hists[j])):
			backgroundhisttemp.append(TH1F("background_"+hists[j][k],"Background fit for "+hists[j][k][hists[j][k].find("2")+1:len(hists[j][k])]+";Mass_{#gamma#gamma} (GeV);Counts",80,80,160))
		backgroundhist.append(backgroundhisttemp)
	print "######################################################################################################################################################################################"
	print "Opening file: %s" %(SignalFiles[i])
	data.append(TFile(pwd+SignalFiles[i]))
	HiggsMass = int(SignalFiles[i][SignalFiles[i].find("Analysis")+8:SignalFiles[i].find("GeV")])
	UpperBin=int(ceil(float(HiggsMass)*1.02)-79)
	LowerBin=int(floor(float(HiggsMass)*0.98)-79)
	#UpperBin=HiggsMass+5-79
	#LowerBin=HiggsMass-5-79
	print "The Mass Range is: %d to %d" %(HiggsMass*0.98,HiggsMass*1.02)
	print "Mass Range is bin number %d to %d" %(LowerBin,UpperBin)
	for j in range(len(BackgroundFiles)):
		#print "Looking at Background File: %s" %BackgroundFiles[j]
		background[j].cd()
		backgroundlist=[]
		for k in range(len(hists)):
			events=[]
			hist=[]
			for l in range(len(hists[k])):
				hist.append(gDirectory.Get(hists[k][l]))
				backgroundhist[k][l].Add(gDirectory.Get(hists[k][l]))
				#import pdb; pdb.set_trace()
				events.append(hist[l].Integral(LowerBin,UpperBin))
				#print "The number of events in %s is: %0.02f" %(hists[k][l],events[l])
			backgroundlist.append(events)
		backgroundmatrix.append(backgroundlist)
	data[i].cd()
	for j in range(len(backgroundsum)):
		for k in range(len(backgroundsum[k])):
			#can = TCanvas("Plots","Plots",1200,1200)
			#backgroundhist[j][k].Draw()
			fit.SetParameters(backgroundhist[j][k].GetMaximum(),-0.02)
			backgroundhist[j][k].Fit(fit,"QLMN","",100,160)
			backgroundintegral[j][k] = fit.Integral(HiggsMass*0.98,HiggsMass*1.02)
			backgroundintegralerror[j][k] = fit.IntegralError(HiggsMass*0.98,HiggsMass*1.02)
			#fit.Draw("same")
			#can.SaveAs(hists[j][k]+".gif")
			for l in range(len(backgroundmatrix)):
				backgroundsum[j][k] = backgroundsum[j][k] + backgroundmatrix[l][j][k]
	for j in range(len(hists)):
		#print "Looking at file %s: %s" %(SignalFiles[i],hists[j][0][hists[j][0].find("2"):hists[j][0].find("All")])
		events=[]
		hist=[]
		for k in range(len(hists[j])):
			hist.append(gDirectory.Get(hists[j][k]))
			events.append(hist[k].Integral(LowerBin,UpperBin))
			#print "The number of events in %s is: %0.02f with %.02f background events so the purity is %.02f%%" %(hists[j][k],events[k],backgroundsum[j][k],100*events[k]/(events[k]+backgroundsum[j][k]))
		eventmatrix.append(events)
	print "\nResults Table for %s" %SignalFiles[i]
	print "Selection\t\t%s\t%s\t\t%s\t\t%s\t\t\t%s" %(SignalFiles[i],BackgroundFiles[0],BackgroundFiles[1],BackgroundFiles[2],BackgroundFiles[3])
	for j in range(len(eventmatrix)):
		if (len(hists[j][0][8:hists[j][0].find("All")])<6):
			print "%s:\t\t\t%.03e (%.03e)\t\t%.03e (%.03e)\t\t%.03e (%.03e)\t\t\t%.03e (%.03e)\t\t%.03e (%.03e)" %(hists[j][0][8:hists[j][0].find("All")],eventmatrix[j][2],eventmatrix[j][3],backgroundmatrix[0][j][2],backgroundmatrix[0][j][3],backgroundmatrix[1][j][2],backgroundmatrix[1][j][3],backgroundmatrix[2][j][2],backgroundmatrix[2][j][3],backgroundmatrix[3][j][2],backgroundmatrix[3][j][3])
		else:
			print "%s:\t\t%.03e (%.03e)\t\t%.03e (%.03e)\t\t%.03e (%.03e)\t\t\t%.03e (%.03e)\t\t%.03e (%.03e)" %(hists[j][0][8:hists[j][0].find("All")],eventmatrix[j][2],eventmatrix[j][3],backgroundmatrix[0][j][2],backgroundmatrix[0][j][3],backgroundmatrix[1][j][2],backgroundmatrix[1][j][3],backgroundmatrix[2][j][2],backgroundmatrix[2][j][3],backgroundmatrix[3][j][2],backgroundmatrix[3][j][3])
	print "Likelihood Sum: " 
	print "\nPurity Results for %s" %SignalFiles[i]
	print "Purity Table:\t\tSignal Events\t\t\tIntegrated Background\t\tIntegrated Error\t\tPurity\t\t\tLog Likelihood\t\tSensitivity\t\tEfficiency X Acceptance"
	for j in range(len(EventReference)):
		if EventReference[i][0]==SignalFiles[i]: NumEvents = EventReference[i][1]
	likelihoodsum = 0
	for j in range(len(eventmatrix)):
		histname = hists[j][0][8:hists[j][0].find("All")]
		purityEB = 100*eventmatrix[j][2]/(backgroundintegral[j][2]+eventmatrix[j][2])
		purityEE = 100*eventmatrix[j][3]/(backgroundintegral[j][3]+eventmatrix[j][3])
		likelihoodEB = eventmatrix[j][2]-(eventmatrix[j][2]+backgroundintegral[j][2])*log(1+eventmatrix[j][2]/backgroundintegral[j][2])
		likelihoodEE = eventmatrix[j][3]-(eventmatrix[j][3]+backgroundintegral[j][3])*log(1+eventmatrix[j][3]/backgroundintegral[j][3])
		sensitivityEB = sqrt(-2*likelihoodEB)
		sensitivityEE = sqrt(-2*likelihoodEE)
		acceptanceEB = eventmatrix[j][2]/NumEvents*100
		acceptanceEE = eventmatrix[j][3]/NumEvents*100
		if (len(histname)<6):
			print "%s:\t\t\t%.03e (%.03e)\t\t%.03e (%.03e)\t\t%.03e (%.03e)\t\t%.02f%% (%.02f%%)\t\t%.03f (%.03f)\t\t%.03f (%.03f)\t\t%.03f (%.03f)" %(histname,eventmatrix[j][2],eventmatrix[j][3],backgroundintegral[j][2],backgroundintegral[j][3],backgroundintegralerror[j][2],backgroundintegralerror[j][3],purityEB,purityEE,likelihoodEB,likelihoodEE,sensitivityEB,sensitivityEE,acceptanceEB,acceptanceEE)
		else:
			print "%s:\t\t%.03e (%.03e)\t\t%.03e (%.03e)\t\t%.03e (%.03e)\t\t%.02f%% (%.02f%%)\t\t%.03f (%.03f)\t\t%.03f (%.03f)\t\t%.03f (%.03f)" %(histname,eventmatrix[j][2],eventmatrix[j][3],backgroundintegral[j][2],backgroundintegral[j][3],backgroundintegralerror[j][2],backgroundintegralerror[j][3],purityEB,purityEE,likelihoodEB,likelihoodEE,sensitivityEB,sensitivityEE,acceptanceEB,acceptanceEE)
		if (j!=0): likelihoodsum = likelihoodsum + likelihoodEB + likelihoodEE
	print "\nLikelihood of (gammaGolden, gamma1goodconv, gamma1poorconv, gamma2conv, gammaleftover) sum is: %.03f and the sensitivity is: %.03f" %(likelihoodsum,sqrt(-2*likelihoodsum)) 
	print "\n"

print "Done!"
#import pdb; pdb.set_trace()
