#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>


using namespace std;

double f(double u) {
  return 2*u;
}


double f_derivace(double u) {
  return 2.0;
}


double alfa(double t) {
  return 1.0;
}


double beta(double t) {
  return 0.0;
}


double casovy_krok(vector<double>& u, double h) {
  double mf = 0;
  for (int i=0; i<u.size(); i++) {
    mf = max(mf, fabs(f_derivace(u[i])));
  }
  return h/mf;
}


double numericky_tok(double ul, double ur) {
  double a = f_derivace( (ul+ur)/2 );
  if (a>0) {
    return f(ul);
  } else {
    return f(ur);
  }
}


void vypocti_uNove(vector<double>& u, double time, double tau, double h, vector<double>& uNove) {
  int m = u.size();
  vector<double> F(m+1);

  F[0] = f(alfa(time));   // Zde by melo byt rozhodovani dle znamenka f_derivace
  F[m] = f(beta(time));   // - || -

  for (int i=1; i<m; i++) {
    F[i] = numericky_tok(u[i-1], u[i]);
  }

  for (int i=0; i<m; i++) {
    uNove[i] = u[i] - tau/h * ( F[i+1] - F[i] );
  }
  
}


void nastav_pocatecni_podminku(vector<double>& u0, double L) {
  int m = u0.size();
  double h = L / m;
  
  for (int i=0; i<m; i++) {
    double x = h*(i+1/2.0);

    if (x>0.25 && x<0.5) {
      u0[i] = 1.0;
    } else {
      u0[i] = 0.0;
    }
  }
  
}

void uloz_data(vector<double>&u, double L, string jmeno) {
  ofstream vystup(jmeno);   // Starsi standard chce vystup(jmeno.c_str()) 
  int m = u.size();
  double h = L / m;
  
  for (int i=0; i<m; i++) {
    double x = h*(i+1/2.0);

    vystup << x << " " << u[i] << endl;
  }  
}


int main() {

  // Priprava dat
  double L = 1.0;      // Delka intervalu
  int    m = 100;      // Pocet bunek
  double tEnd = 0.1;
  double h = L / m;
  
  vector<double> u(m);
  vector<double> uNove(m);
  
  // Vypocet
  nastav_pocatecni_podminku(u, L);
  double time = 0;
  while (time < tEnd) {
    double tau = 0.8 * casovy_krok(u, h);
    vypocti_uNove(u, time, tau, h, uNove);
    time += tau;
    for (int i=0; i<m; i++) u[i] = uNove[i];
  }
  
  
  // Ulozeni
  uloz_data(u, L, "vysledek.txt");
	    
  return 0;
}

