//          Copyright Jiri Furst 2016
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef EULER_HPP
#define EULER_HPP

#include <cmath>

#include "Vector.hpp"

namespace Euler {

  const size_t MODEL_SIZE = 3;   // Pocet rovnic modelu

  using Vars = Vector<MODEL_SIZE,double>;

  enum { RHO, RHO_U, RHO_E };    // Pojmenovani konzervativnich promennych

  double kappa = 1.4;            // Poissonova konstanta


  inline double density(const Vars& w) { return w[RHO]; }

  inline double velocity(const Vars& w) { return w[RHO_U] / w[RHO]; }

  inline double total_energy(const Vars& w) { return w[RHO_E] / w[RHO]; }


  inline double pressure(const Vars& w) {
    auto rho = density(w);
    auto u   = velocity(w);
    auto e   = total_energy(w);

    return (kappa - 1.0) * rho * ( e - 0.5 * u * u);
  }


  inline double sound_speed(const Vars& w) {
    auto rho = density(w);
    auto p   = pressure(w);

    return sqrt( kappa * p / rho );
  }

  
  inline Vector<MODEL_SIZE,double> wave_speeds(const Vars& w) {
    auto u = velocity(w);
    auto a = sound_speed(w);

    return Vector<MODEL_SIZE,double>( { u-a, u, u+a} );
  }


  inline double spectral_radius(const Vars& w) {
    auto u = velocity(w);
    auto a = sound_speed(w);

    return fabs(u)+a;
  }


  inline Vars flux(const Vars& w) {
    auto u = velocity(w);
    auto p = pressure(w);

    return u * w + Vars({ 0., p, p*u });
  }

};

#endif
