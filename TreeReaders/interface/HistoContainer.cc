#include "ND_Hto2Photons/TreeReaders/interface/HistoContainer.h"
#include <utility>
#include <iostream>

HistoContainer::HistoContainer() 
{}

HistoContainer::~HistoContainer() 
{}

void HistoContainer::Add(string name, int bins, float xmin, float xmax) {
  TH1F temp(name.c_str(), name.c_str(), bins, xmin, xmax);
  temp.Sumw2();
  h1[name] = temp;
}

void HistoContainer::Add(string name, string title, int bins, float xmin, float xmax) {
  TH1F temp(name.c_str(), title.c_str(), bins, xmin, xmax);
  temp.Sumw2();
  h1[name] = temp;
}

void HistoContainer::Add(string name, int binsx, float xmin, float xmax,
			 int binsy, float ymin, float ymax) {
		       
  TH2F temp(name.c_str(), name.c_str(), binsx, xmin, xmax, binsy, ymin, ymax);
  temp.Sumw2();
  h2[name] = temp;
}

void HistoContainer::Add(string name, string title, int binsx, float xmin, float xmax,
			 int binsy, float ymin, float ymax) {
		       
  TH2F temp(name.c_str(), title.c_str(), binsx, xmin, xmax, binsy, ymin, ymax);
  temp.Sumw2();
  h2[name] = temp;
}

void HistoContainer::Add(string name, int binsx, float xmin, float xmax,
			 float ymin, float ymax) {
		       
  TProfile temp(name.c_str(), name.c_str(), binsx, xmin, xmax, ymin, ymax);
  temp.Sumw2();
  hp[name] = temp;
} 

void HistoContainer::Add(string name, string title, int binsx, float xmin, float xmax,
			 float ymin, float ymax) {
		       
  TProfile temp(name.c_str(), title.c_str(), binsx, xmin, xmax, ymin, ymax);
  temp.Sumw2();
  hp[name] = temp;
}

void HistoContainer::Fill(string name, float value) {

  std::map<string, TH1F>::const_iterator it = h1.find(name);
  if (it != h1.end()) {
    h1[name].Fill(value);
    return;
  }

  std::cerr << "ERROR !: histogram " << name << " is not a TH1F." << std::endl;
}

void HistoContainer::Fill(string name, float valuex, float valuey) { 

  std::map<string, TH1F>::const_iterator it = h1.find(name);
  if (it != h1.end()) {
    h1[name].Fill(valuex,valuey);
    return;
  }

  std::map<string, TH2F>::const_iterator it2 = h2.find(name);
  if (it2 != h2.end()) {
    h2[name].Fill(valuex, valuey);
    return;
  }

  std::map<string, TProfile>::const_iterator itp = hp.find(name);
  if (itp != hp.end()) {
    hp[name].Fill(valuex, valuey);
    return;
  }
  
  std::cerr << "ERROR !: histogram " << name << " is nor a TH2F nor a TProfile." << std::endl;
}

void HistoContainer::Fill(string name, float valuex, float valuey, float weight) { 

  std::map<string, TH2F>::const_iterator it2 = h2.find(name);
  if (it2 != h2.end()) {
    h2[name].Fill(valuex, valuey, weight);
    return;
  }

  std::map<string, TProfile>::const_iterator itp = hp.find(name);
  if (itp != hp.end()) {
    hp[name].Fill(valuex, valuey, weight);
    return;
  }
  
  std::cerr << "ERROR !: histogram " << name << " is nor a TH2F" << std::endl;
}

void HistoContainer::Fill(string name, string addon, float value) {

  name.append(addon);
  std::map<string, TH1F>::const_iterator it = h1.find(name);
  if (it != h1.end()) {
    h1[name].Fill(value);
    return;
  }

  std::cerr << "ERROR !: histogram " << name << " is not a TH1F." << std::endl;
}

void HistoContainer::Fill(string name, string addon, float valuex, float valuey) { 

  name.append(addon);
  std::map<string, TH1F>::const_iterator it = h1.find(name);
  if (it != h1.end()) {
    h1[name].Fill(valuex,valuey);
    return;
  }
  
  std::map<string, TH2F>::const_iterator it2 = h2.find(name);
  if (it2 != h2.end()) {
    h2[name].Fill(valuex, valuey);
    return;
  }

  std::map<string, TProfile>::const_iterator itp = hp.find(name);
  if (itp != hp.end()) {
    hp[name].Fill(valuex, valuey);
    return;
  }
  
  std::cerr << "ERROR !: histogram " << name << " is nor a TH2F nor a TProfile." << std::endl;
}

void HistoContainer::Fill(string name, string addon, float valuex, float valuey, float weight) { 

  name.append(addon);
  std::map<string, TH2F>::const_iterator it2 = h2.find(name);
  if (it2 != h2.end()) {
    h2[name].Fill(valuex, valuey, weight);
    return;
  }

  std::map<string, TProfile>::const_iterator itp = hp.find(name);
  if (itp != hp.end()) {
    hp[name].Fill(valuex, valuey, weight);
    return;
  }

  std::cerr << "ERROR !: histogram " << name << " is nor a TH2F" << std::endl;
}


void HistoContainer::Fill(string name, string addon1, string addon2, float value) {

  name.append(addon1);
  name.append(addon2);
  std::map<string, TH1F>::const_iterator it = h1.find(name);
  if (it != h1.end()) {
    h1[name].Fill(value);
    return;
  }

  std::cerr << "ERROR !: histogram " << name << " is not a TH1F." << std::endl;
}

void HistoContainer::Fill(string name, string addon1, string addon2, float valuex, float valuey) { 

  name.append(addon1);
  name.append(addon2);
  std::map<string, TH1F>::const_iterator it = h1.find(name);
  if (it != h1.end()) {
    h1[name].Fill(valuex,valuey);
    return;
  }
  
  std::map<string, TH2F>::const_iterator it2 = h2.find(name);
  if (it2 != h2.end()) {
    h2[name].Fill(valuex, valuey);
    return;
  }

  std::map<string, TProfile>::const_iterator itp = hp.find(name);
  if (itp != hp.end()) {
    hp[name].Fill(valuex, valuey);
    return;
  }
  
  std::cerr << "ERROR !: histogram " << name << " is nor a TH2F nor a TProfile." << std::endl;
}

void HistoContainer::Fill(string name, string addon1, string addon2, float valuex, float valuey, float weight) { 

  name.append(addon1);
  name.append(addon2);
  std::map<string, TH2F>::const_iterator it2 = h2.find(name);
  if (it2 != h2.end()) {
    h2[name].Fill(valuex, valuey, weight);
    return;
  }

  std::map<string, TProfile>::const_iterator itp = hp.find(name);
  if (itp != hp.end()) {
    hp[name].Fill(valuex, valuey, weight);
    return;
  }

  std::cerr << "ERROR !: histogram " << name << " is nor a TH2F" << std::endl;
}



double HistoContainer::UpperLimit(string name) {

  std::map<string, TH1F>::const_iterator it = h1.find(name);
  if (it != h1.end()) {
    return h1[name].GetBinLowEdge(h1[name].GetNbinsX()+1);
  }

  std::cerr << "ERROR !: histogram " << name << " is not a TH1F." << std::endl;
  return -999;

}

double HistoContainer::LowerLimit(string name) {

  std::map<string, TH1F>::const_iterator it = h1.find(name);
  if (it != h1.end()) {
    return h1[name].GetBinLowEdge(1);
  }

  std::cerr << "ERROR !: histogram " << name << " is not a TH1F." << std::endl;
  return -999;
}

void HistoContainer::Save() {
  std::map<string, TH1F>::const_iterator it;
  for (it = h1.begin(); it != h1.end(); ++it)
    it->second.Write();

  std::map<string, TH2F>::const_iterator it2;
  for (it2 = h2.begin(); it2 != h2.end(); ++it2)
    it2->second.Write(); 

  std::map<string, TProfile>::const_iterator itp;
  for (itp = hp.begin(); itp != hp.end(); ++itp)
    itp->second.Write();
}
