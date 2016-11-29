//          Copyright Jiri Furst 2016
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "Vector.hpp"

template <size_t N, size_t M = N, typename T=double>
class Matrix {
public:
  typedef Vector<M,T> Row;

  Matrix() = default;
  Matrix(const Matrix& x) = default;

  inline size_t size() const { return N; }

  Row& operator[](size_t i) { return data_[i]; }
  const Row& operator[](size_t i) const { return data_[i]; }
 
private:
  Vector<N,Row> data_;;

};


// a + b
template <size_t N, size_t M, typename T>
inline Matrix<N,M,T> operator+(const Matrix<N,M,T>& u, const Matrix<N,M,T>& v) {
  Matrix<N,M,T> res;
  for (size_t i=0; i<N; i++)
    res[i] = u[i] + v[i];
  return res;
}

// a - b
template <size_t N, size_t M, typename T>
inline Matrix<N,M,T> operator-(const Matrix<N,M,T>& u, const Matrix<N,M,T>& v) {
  Matrix<N,M,T> res;
  for (size_t i=0; i<N; i++)
    res[i] = u[i] - v[i];
  return res;
}

// a * u
template <size_t N, size_t M, typename T>
inline Matrix<N,M,T> operator*(T a, const Matrix<N,M,T>& u) {
  Matrix<N,M,T> res;
  for (size_t i=0; i<N; i++)
    res[i] = a * u[i];
  return res;
}

// u / a
template <size_t N, size_t M, typename T>
inline Matrix<N,M,T> operator/(const Matrix<N,M,T>& u, T a) {
  Matrix<N,M,T> res;
  for (size_t i=0; i<N; i++)
    res[i] = u[i] / a;
  return res;
}


// a * x
template <size_t N, size_t M, typename T>
inline Vector<N,T> operator*(const Matrix<N,M,T>& a, const Vector<M,T>& x) {
  Vector<N,T> res;
  for (size_t i=0; i<N; i++)
    res[i] = dot(a[i], x);
  return res;
}

#endif
