//          Copyright Jiri Furst 2017
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <fstream>
#include <vector>
#include <limits>

#include "Euler.hpp"

using namespace std;
using namespace Euler;


double area(double x) {
  return 1 + pow(x-0.5, 2);
}


Vars HLL(const Vars& wl, const Vars& wr) {

  auto lambdaL = wave_speeds(wl);
  auto lambdaR = wave_speeds(wr);

  auto sl = min(lambdaL[0],lambdaR[0]);
  auto sr = max(lambdaL[2],lambdaR[2]);

  if ( 0 <= sl )
    return flux(wl);
  else if ( 0 < sr )
    return ( sr * flux(wl) - sl * flux(wr) + sl*sr * (wr - wl) ) / (sr - sl);
  else
    return flux(wr);
}


void save_data(const char* filename, const vector<Vars>& w, double dx) {
  ofstream f(filename);
  auto x = dx/2;
  for (auto wi : w) {
    f << x;
    for (auto j=0; j<wi.size(); j++) f << "\t" << wi[j];
    f << "\t" << pressure(wi) << "\t" << velocity(wi)/sound_speed(wi);
    f << endl;
    x += dx;
  }
}


double maximal_timestep(const vector<Vars>& w, double dx) {
  double dt = numeric_limits<double>::max();
  for (auto wi : w) 
    dt = min( dt, dx / spectral_radius(wi) );
  return dt;
}


int main() {

  double p0 = 1.e5;
  double T0 = 293.15;
  double r  = 287.05;
  double p2 = 0.75*p0;
  
  const int n = 400;
  double dx = 1.0 / n;

  
  vector<Vars> w(n);
  for (auto i=0; i<w.size(); i++)
    w[i] = Vars({ p0 / (r*T0), 0.0, p0/(kappa-1) });
  
  save_data("w0.dat", w, dx);

  vector<Vars> wn(w.size());
  vector<Vars> f(w.size()+1);

  for (auto iter = 0; iter < 20000; iter++) {

    // Okrajove podminky na vstupu
    auto pIn   = min(pressure(w[0]), p0);
    auto M2    = (pow(pIn/p0, (1-kappa)/kappa) - 1) * 2/(kappa-1);
    auto rhoIn = p0/(r*T0) * pow(1 + (kappa-1)/2*M2, 1/(1-kappa));
    auto uIn   = sqrt(M2 * kappa*pIn/rhoIn);
    Vars Win   = Vars({ rhoIn, rhoIn*uIn, pIn/(kappa-1) + 0.5*rhoIn*uIn*uIn });
    f[0] = flux(Win);

    // Okrajove podminky na vystupu
    Vars Wout = w[n-1];
    Wout[RHO_E] = p2/(kappa-1) + 0.5*pow(Wout[RHO_U],2)/Wout[RHO];
    f[n] = flux(Wout);

    for (auto i=1; i<w.size(); i++)
      f[i] = HLL(w[i-1], w[i]);


    auto dt = 0.4 * maximal_timestep(w, dx);
    
    for (auto i=0; i<w.size(); i++) {
      auto Al = area(i*dx);
      auto Ar = area((i+1)*dx);
      auto Ac = (Al + Ar) / 2;
      
      Vars Q = Vars({ 0.0, (Ar - Al)/dx * pressure(w[i]), 0.0});
      
      wn[i] = w[i] - dt/(dx*Ac) * (Ar*f[i+1] - Al*f[i]) + dt/Ac * Q;
    }
    
    for (auto i=0; i<w.size(); i++) w[i] = wn[i];

    if (iter % 100 == 0)
      cout << iter << "\t" << dt << endl;
  }
    
  save_data("w.dat", w, dx);

  return 0;
}
