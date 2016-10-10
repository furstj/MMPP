// Copyright 2016 Jiri Furst
// Distributed under MIT license
// See file ../LICENSE for details

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cassert>

using namespace std;

//- Initial condition
double u0(double x) 
{
  if (x<0.25 || x>0.5) return 0.0;
  else return cos( (8*x-1)*M_PI) + 1;
}

//- Save result to a file
void save_result(const string& filename, const vector<double>& x, const vector<double> u)
{
  assert( x.size() == u.size() );

  ofstream f(filename.c_str());
  for (size_t i=0; i<x.size(); i++)
    f << x[i] << "\t" << u[i] << endl;
}

int main()
{
  const int   nx = 101;         // Number of points 
  const double L = 1.0;         // Length of the domain
  const double a = 1.0;         // Advection speed
  const double Tend = 0.25;

  vector<double> x(nx);
  vector<double> u(nx);
  
  for (size_t i=0; i<x.size(); i++) {
    x[i] = i * double(L)/(nx-1);
    u[i] = u0(x[i]);
  };

  save_result("u0.dat", x, u);

  double dx = L / (nx-1);
  double dt = 0.8 * dx / fabs(a);

  int niter = ceil(Tend/dt);
  dt = Tend / niter;

  vector<double> un(nx);

  double t = 0;
  for (int n=0; n<niter; n++) {

    for (int i=1; i<nx; i++) 
      un[i] = u[i] - a*dt/dx * (u[i] - u[i-1]);

    t += dt;
    for (int i=1; i<nx; i++) u[i] = un[i];
  }

  save_result("u.dat", x, u);
  
  return 0;
}

