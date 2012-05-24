#ifndef HISTOCONTAINER
#define HISTOCONTAINER

#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <map>
#include <string>

class HistoContainer {

public:
  HistoContainer();
  ~HistoContainer();

  void Add(string name, int bins, float xmin, float xmax);
  void Add(string name, string title, int bins, float xmin, float xmax);
  void Add(string name, int binsx, float xmin, float xmax, int binsy, float ymin, float ymax);
  void Add(string name, string title, int binsx, float xmin, float xmax, int binsy, float ymin, float ymax);
  void Add(string name, int binsx, float xmin, float xmax, float ymin, float ymax);
  void Add(string name, string title, int binsx, float xmin, float xmax, float ymin, float ymax);

  void Fill(string name, float value, bool overflow);
  void Fill(string name, float valuex, float valuey, bool overflow);
  void Fill(string name, float valuex, float valuey, float weight);

  void Fill(string name, string addon, float value, bool overflow);
  void Fill(string name, string addon, float valuex, float valuey, bool overflow);
  void Fill(string name, string addon, float valuex, float valuey, float weight);

  void Fill(string name, string addon1, string addon2, float value, bool overflow);
  void Fill(string name, string addon1, string addon2, float valuex, float valuey, bool overflow);
  void Fill(string name, string addon1, string addon2, float valuex, float valuey, float weight);

  double UpperLimit(string name);
  double LowerLimit(string name);

  void Save();
  
 private:
  std::map<string, TH1F> h1;
  std::map<string, TH2F> h2;
  std::map<string, TProfile> hp;
};

#endif
