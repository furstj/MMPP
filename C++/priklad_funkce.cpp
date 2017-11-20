#include <iostream>         // Pouzivej knihovnu pro vstupy a vystupy (klavesnice, obrazovka)
#include <cmath>            // Pouzivej matematickou knihovnu

using namespace std;        // Hledej objekty v knihovne std (std::vector -> vector)

// Vypocet funkce f(x) = sin(sqrt(x+1))

double f(double x) {     // Predani x hodnotou (parametr se kopiruje do x)
  double vysledek;

  x = x + 1;
  vysledek = sin(sqrt(x));

  return vysledek;
}


double g(double& x) {     // Predani x referenci (x zastupuje nejakou promennou)
  double vysledek;

  x = x + 1;
  vysledek = sin(sqrt(x));

  return vysledek;
}


int main(int argc, char** argv) {

  double z = 1.5;

  double y = f(z);

  cout << "Vysledek vypoctu je " << y << endl;
  cout << "z je " << z << endl;


  double v = 1.5;
  double u = g(v);

  cout << "Vysledek vypoctu je " << u << endl;
  cout << "v je " << v << endl;
  
  return 0;
}
