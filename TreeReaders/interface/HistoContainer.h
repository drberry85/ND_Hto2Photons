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
    
  void Add(string, int, float, float);
  void Add(string, string, int, float, float);
  void Add(string, int, float, float, int, float, float);
  void Add(string, string, int, float, float, int, float, float);
  void Add(string, int, float, float, float, float);

  void Fill(string, float);
  void Fill(string, float, float);
  void Fill(string, float, float, float);

  void Fill(string, string, float);
  void Fill(string, string, float, float);
  void Fill(string, string, float, float, float);
  
  void Save();
  
 private:
  std::map<string, TH1F> h1;
  std::map<string, TH2F> h2;
  std::map<string, TProfile> hp;
};

#endif
