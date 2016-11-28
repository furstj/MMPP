//          Copyright Jiri Furst 2016
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <array>

template <size_t N, typename T=double>
class Vector {
public:

  Vector() = default;
  Vector(const Vector& x) = default;
  Vector(const std::array<T,N>& x) : data_(x) {};

  inline size_t size() const { return N; }

  T& operator[](size_t i) { return data_[i]; }
  const T& operator[](size_t i) const { return data_[i]; }
 
private:
  std::array<T,N> data_;

};


// u + v
template <size_t N, typename T>
inline Vector<N,T> operator+(const Vector<N,T>& u, const Vector<N,T>& v) {
  Vector<N,T> res;
  for (size_t i=0; i<N; i++)
    res[i] = u[i] + v[i];
  return res;
}

// u - v
template <size_t N, typename T>
inline Vector<N,T> operator-(const Vector<N,T>& u, const Vector<N,T>& v) {
  Vector<N,T> res;
  for (size_t i=0; i<N; i++)
    res[i] = u[i] - v[i];
  return res;
}

// a * u
template <size_t N, typename T>
inline Vector<N,T> operator*(T a, const Vector<N,T>& u) {
  Vector<N,T> res;
  for (size_t i=0; i<N; i++)
    res[i] = a * u[i];
  return res;
}

// u / a
template <size_t N, typename T>
inline Vector<N,T> operator/(const Vector<N,T>& u, T a) {
  Vector<N,T> res;
  for (size_t i=0; i<N; i++)
    res[i] = u[i] / a;
  return res;
}

template <size_t N, typename T>
inline T dot(const Vector<N,T>& u, const Vector<N,T>&v) {
  T res = 0;
  for (size_t i=0; i<N; i++)
    res += u[i] * v[i];
  return res;
}


