#include <iostream>         // Pouzivej knihovnu pro vstupy a vystupy (klavesnice, obrazovka)
#include <vector>
#include <string>

#include "prace_se_soubory.hpp"


using namespace std;        // Hledej objekty v knihovne std (std::vector -> vector)


int main(int argc, char** argv) {

  vector<double> p = nacti_pole("vstup.txt");  // Soubor vstup.txt obsahuje delku a prvky pole

  for (int i=0; i<p.size(); i++) {
    p[i] = 2*p[i] + 1;
  }

  uloz_pole(p, "vystup.txt");
  
  return 0;
}
