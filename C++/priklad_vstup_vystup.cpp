#include <vector>           // Pouzivej knihovnu pro pole
#include <iostream>         // Pouzivej knihovnu pro vstupy a vystupy (klavesnice, obrazovka)

using namespace std;        // Hledej objekty v knihovne std (std::vector -> vector)


// Nacteme pocet prvku z klavesnice, vypocteme prvky pole, vytiskneme

int main(int argc, char** argv) {

  cout << "Zadej delku pole" << endl;

  int n;
  cin >> n;
  cout << "Zadana delka je " << n << endl;
  
  vector<double> p(n);
  
  for (int i=0; i<p.size(); i++) {   // p.size() vraci delku pole
    p[i] = 0.5*i + 1;
  }

  for (int i=0; i<p.size(); i++) {
    cout << "p[" << i << "] = " << p[i] << endl;
  }
  
  return 0;
}
