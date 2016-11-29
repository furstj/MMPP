//          Copyright Jiri Furst 2016
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


Vars initial_condition(double x) {
  if (x < 0.5) 
    return Vars( {10, 0, 1.e6/kappa} );
  else
    return Vars( {1, 0, 1.e5/kappa} );
}


void save_data(const char* filename, const vector<Vars>& w, double dx) {
  ofstream f(filename);
  auto x = dx/2;
  for (auto wi : w) {
    f << x;
    for (auto j=0; j<wi.size(); j++) f << "\t" << wi[j];
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
  const int n = 400;
  double dx = 1.0 / n;
  
  vector<Vars> w(n);
  for (auto i=0; i<w.size(); i++)
    w[i] = initial_condition( (i+0.5)*dx );
  save_data("w0.dat", w, dx);

  vector<Vars> wn(w.size());
  vector<Vars> f(w.size()+1);

  for (auto iter = 0; iter < n; iter++) {

    // Osizene okrajove podminky
    f[0] = HLL(w[0], w[0]);
    f[n] = HLL(w[n-1], w[n-1]);

    for (auto i=1; i<w.size(); i++)
      f[i] = HLL(w[i-1], w[i]);


    auto dt = 0.4 * maximal_timestep(w, dx);
    for (auto i=0; i<w.size(); i++)
      wn[i] = w[i] - dt/dx * (f[i+1] - f[i]);

    for (auto i=0; i<w.size(); i++) w[i] = wn[i];

  }
    
  save_data("w.dat", w, dx);
    
  return 0;
}
