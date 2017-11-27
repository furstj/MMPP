#include <iostream>
#include <vector>


using namespace std;


struct Dvojice {
  double eta;
  double q;
};


Dvojice operator+(Dvojice a, Dvojice b) {
  Dvojice soucet;
  soucet.eta = a.eta + b.eta;
  soucet.q   = a.q   + b.q;
  return soucet;
}


Dvojice operator-(Dvojice a, Dvojice b) {
  Dvojice rozdil;
  rozdil.eta = a.eta - b.eta;
  rozdil.q   = a.q   - b.q;
  return rozdil;
}


Dvojice operator*(double a, Dvojice b) {
  Dvojice vysledek;
  vysledek.eta = a * b.eta;
  vysledek.q   = a * b.q;
  return vysledek;
}


int main() {

  Dvojice x;
  x.eta = 1;
  x.q   = 2;
  

  vector<Dvojice> u(100);
  for (int i=0; i<u.size(); i++) {
    u[i].eta = 1.;
    u[i].q   = 2.;
  }

  Dvojice y, z;
  y = x;

  /*
  z.eta = x.eta + y.eta;
  z.q   = x.q   + y.q;
  */

  z = x + y;

  z = x - y;

  Dvojice p = 2.0 * x;

  return 0;
}
