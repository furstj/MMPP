#include <vector>           // Pouzivej knihovnu pro pole
#include <iostream>         // Pouzivej knihovnu pro vstupy a vystupy (klavesnice, obrazovka)
#include <fstream>          // Pouzivej knihovnu pro vstupy a vystupy (soubor)

using namespace std;        // Hledej objekty v knihovne std (std::vector -> vector)



int main(int argc, char** argv) {

  int n = 10;
  vector<double> p(n);
  
  for (int i=0; i<p.size(); i++) {   // p.size() vraci delku pole
    p[i] = 0.5*i + 1;
  }

  ofstream vystup("pole.dat");       // Vytvor soubor "pole.dat" (pro zapis)
  
  for (int i=0; i<p.size(); i++) {
    vystup << p[i] << endl;
  }

  vystup.close();


  // Nacteni dat do pole q
  vector<double> q(n);
  ifstream vstup("pole.dat");         // Otevri pole.dat ke cteni

  for (int i=0; i<q.size(); i++) {
    vstup >> q[i];
  }

  // Tisk na obrazovku
  cout << "Pole q" << endl;
  for (int i=0; i<q.size(); i++) {
      cout << q[i] << endl;
  }
  
  return 0;
}
