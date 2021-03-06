#pragma once

/** @file Point.hpp
 * @brief Define the Point class for 3D points.
 */

#include <iostream>
#include <cmath>
#include "../Util/TreeTypes.hpp"

#define for_i for(std::size_t i = 0; i != 3; ++i)

/** @class Point
 * @brief Class representing 3D points and vectors.
 *
 * Point contains methods that support use as points in 3D space, and
 * use as 3-dimensional vectors (that can, for example, be cross-producted).
 *
 * Its x, y, and z components are publicly accessible under those names. They
 * can also be accessed as coordinate[0], coordinate[1], and coordinate[2],
 * respectively.
 */
struct Point {
  union {
    struct {
      real_type x;
      real_type y;
      real_type z;
    };
    real_type elem[3];
  };

  typedef real_type          value_type;
  typedef real_type&         reference;
  typedef const real_type&   const_reference;
  typedef real_type*         iterator;
  typedef const real_type*   const_iterator;
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;

  // CONSTRUCTORS

  Point() {
    for_i elem[i] = value_type();
  }
  explicit Point(value_type b) {
    for_i elem[i] = b;
  }
  Point(value_type b0, value_type b1, value_type b2) {
    elem[0] = b0; elem[1] = b1; elem[2] = b2;
  }

  // Ad hoc constructors
  Point(const FloatArray<3>& b) {
    for_i elem[i] = b[i];
  }
  Point(const DoubleArray<3>& b) {
    for_i elem[i] = b[i];
  }

  // COMPARATORS

  bool operator==(const Point& b) const {
    for_i if (elem[i] != b[i]) return false;
    return true;
  }
  bool operator!=(const Point& b) const {
    return !(*this == b);
  }

  // MODIFIERS

  /** Add scalar @a b to this Point */
  Point& operator+=(value_type b) {
    for_i elem[i] += b;
    return *this;
  }
  /** Subtract scalar @a b from this Point */
  Point& operator-=(value_type b) {
    for_i elem[i] -= b;
    return *this;
  }
  /** Scale this Point up by scalar @a b */
  Point& operator*=(value_type b) {
    for_i elem[i] *= b;
    return *this;
  }
  /** Scale this Point down by scalar @a b */
  Point& operator/=(value_type b) {
    for_i elem[i] /= b;
    return *this;
  }
  /** Add Point @a b to this Point */
  Point& operator+=(const Point& b) {
    for_i elem[i] += b[i];
    return *this;
  }
  /** Subtract Point @a b from this Point */
  Point& operator-=(const Point& b) {
    for_i elem[i] -= b[i];
    return *this;
  }
  /** Scale this Point up by factors in @a b */
  Point& operator*=(const Point& b) {
    for_i elem[i] *= b[i];
    return *this;
  }
  /** Scale this Point down by factors in @a b */
  Point& operator/=(const Point& b) {
    for_i elem[i] /= b[i];
    return *this;
  }

  // ACCESSORS

  reference       operator[](size_type i)       { return elem[i]; }
  const_reference operator[](size_type i) const { return elem[i]; }

  iterator        data()       { return elem; }
  const_iterator  data() const { return elem; }

  reference       front()       { return elem[0]; }
  const_reference front() const { return elem[0]; }
  reference       back()        { return elem[2]; }
  const_reference back()  const { return elem[2]; }

  static constexpr size_type     size() { return 3; }
  static constexpr size_type max_size() { return 3; }
  static constexpr bool         empty() { return false; }

  // ITERATORS

  iterator        begin()       { return elem; }
  const_iterator  begin() const { return elem; }
  const_iterator cbegin() const { return elem; }

  iterator          end()       { return elem+3; }
  const_iterator    end() const { return elem+3; }
  const_iterator   cend() const { return elem+3; }
};

// STREAM OPERATORS

/** Write a Point to an output stream */
std::ostream& operator<<(std::ostream& s, const Point& a) {
  return (s << a.x << ' ' << a.y << ' ' << a.z);
}
/** Read a Vec from an input stream */
std::istream& operator>>(std::istream& s, Point& a) {
  return (s >> a.x >> a.y >> a.z);
}

// ARITHMETIC OPERATORS

/** Unary negation: Return -@a a */
Point operator-(const Point& a) {
  return Point(-a.x, -a.y, -a.z);
}
/** Unary plus: Return @a a. ("+a" should work if "-a" works.) */
Point operator+(const Point& a) {
  return a;
}
Point operator+(Point a, const Point& b) {
  return a += b;
}
Point operator+(Point a, real_type b) {
  return a += b;
}
Point operator+(real_type b, Point a) {
  return a += b;
}
Point operator-(Point a, const Point& b) {
  return a -= b;
}
Point operator-(Point a, real_type b) {
  return a -= b;
}
Point operator-(real_type b, const Point& a) {
  return (-a) += b;
}
Point operator*(Point a, const Point& b) {
  return a *= b;
}
Point operator*(Point a, real_type b) {
  return a *= b;
}
Point operator*(real_type b, Point a) {
  return a *= b;
}
Point operator/(Point a, const Point& b) {
  return a /= b;
}
Point operator/(Point a, real_type b) {
  return a /= b;
}

// NORMS AND MATH OPERATORS

/** Compute cross product of two 3D Points */
Point cross(const Point& a, const Point& b) {
  return Point(a[1]*b[2] - a[2]*b[1],
               a[2]*b[0] - a[0]*b[2],
               a[0]*b[1] - a[1]*b[0]);
}

/** Compute the inner product of two Points */
real_type inner_prod(const Point& a, const Point& b) {
  real_type v = 0;
  for_i v += a[i]*b[i];
  return v;
}
real_type dot(const Point& a, const Point& b) {
  return inner_prod(a,b);
}

/** Compute the squared L2 norm of a Point */
real_type normSq(const Point& a) {
  real_type v = 0;
  for_i v += a[i]*a[i];
  return v;
}
/** Compute the L2 norm of a Point */
real_type norm(const Point& a) {
  return std::sqrt(normSq(a));
}
real_type norm_2(const Point& a) {
  return norm(a);
}
/** Compute the L1 norm of a Point */
real_type norm_1(const Point& a) {
  real_type v = 0;
  for_i v += std::abs(a[i]);
  return v;
}
/** Compute the L-infinity norm of a Point */
real_type norm_inf(const Point& a) {
  real_type v = 0;
  for_i v = std::max(v, std::abs(a[i]));
  return v;
}

#undef for_i
