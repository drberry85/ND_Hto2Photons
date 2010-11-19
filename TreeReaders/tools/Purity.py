print "Loading Root..."

from ROOT import *
gROOT.Macro("$HOME/rootlogon.C")
gStyle.SetOptStat(000000)

print "Setting Initial Parameters."
pwd = "/data/ndpc2/c/HiggsGammaGamma/CMSSW_3_8_5_patch3/src/ND_Hto2Photons/TreeReaders/"
SignalFiles = ["HiggsAnalysis110GeV.root","HiggsAnalysis120GeV.root","HiggsAnalysis150GeV.root"]
BackgroundFiles = ["PhotonPlusJet.root","Born.root","Box.root"]
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
	print "Opening file: %s" %(SignalFiles[i])
	data.append(TFile(pwd+SignalFiles[i]))
	if SignalFiles[i].find("90GeV")!=-1 or SignalFiles[i].find("110GeV")!=-1 or SignalFiles[i].find("120GeV")!=-1 or SignalFiles[i].find("150GeV")!=-1:
		if SignalFiles[i].find("90GeV")!=-1: HiggsMass=90;
		if SignalFiles[i].find("110GeV")!=-1: HiggsMass=110;
		if SignalFiles[i].find("120GeV")!=-1: HiggsMass=120;
		if SignalFiles[i].find("150GeV")!=-1: HiggsMass=150;
		LowerBin=int(floor(HiggsMass*0.98)-80)
		UpperBin=int(floor(HiggsMass*1.02)+1-80)
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
		for j in range(len(backgroundmatrix)):
			if j==0: backgroundsum=backgroundmatrix[j]
			else:
				for k in range(len(backgroundmatrix[j])):
					for l in range(len(backgroundmatrix[j][k])):
						backgroundsum[k][l] = backgroundsum[k][l] + backgroundmatrix[j][k][l]
		for j in range(len(hists)):
			print "Looking at file %s: %s" %(SignalFiles[i],hists[j][0][hists[j][0].find("_mass_"):hists[j][0].find("All")])
			events=[]
			hist=[]
			for k in range(len(hists[j])):
				hist.append(gDirectory.Get(hists[j][k]))
				events.append(hist[k].Integral(LowerBin,UpperBin))
				print "The number of events in %s is: %0.02f with %.02f background events so the purity is %.02f%%" %(hists[j][k],events[k],backgroundsum[j][k],100*events[k]/(events[k]+backgroundsum[j][k]))
			
