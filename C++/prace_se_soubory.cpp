#include "prace_se_soubory.hpp"

#include <vector>
#include <fstream>

using namespace std;


vector<double> nacti_pole(string jmeno) {
  ifstream vstup(jmeno);
  int n;
  vstup >> n;
  vector<double> p(n);
  for (int i=0; i<n; i++) {
    vstup >> p[i];
  }

  return p;
}


void uloz_pole(std::vector<double> x, std::string jmeno) {
  ofstream vystup(jmeno);
  vystup <<  x.size() << endl;
  for (int i=0; i<x.size(); i++) {
    vystup << x[i] << endl;
  }
}
