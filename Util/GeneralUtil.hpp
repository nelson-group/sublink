#pragma once
/** @file GeneralUtil.hpp
 * @brief General utilities, some of them taken from CS207.
 */

#include <iostream>
#include <vector>

////////////////////////////
// Some tricks from CS207 //
////////////////////////////

/** Derive operator>, operator<=, and operator>= from a class's operator<.
 *
 * Usage:
 * class MyClass : private less_than_comparable<MyClass> {
 *   bool operator<(const MyClass& c) { ... }
 * };
 */
template <typename T>
struct less_than_comparable {
  friend bool operator> (const T& a, const T& b) { return   b < a;  }
  friend bool operator<=(const T& a, const T& b) { return !(b < a); }
  friend bool operator>=(const T& a, const T& b) { return !(a < b); }
};

/** Derive operator!= from a class's operator==.
 *
 * Usage:
 * class MyClass : private equality_comparable<MyClass> {
 *   bool operator==(const MyClass& c) { ... }
 * };
 */
template <typename T>
struct equality_comparable {
  friend bool operator!=(const T& a, const T& b) { return !(a == b); }
};

/** Derive operator!=, operator>, operator<=, and operator>=
 * from a class's operator< and operator==.
 *
 * Usage:
 * class MyClass : private equality_comparable<MyClass> {
 *   bool operator< (const MyClass& c) { ... }
 *   bool operator==(const MyClass& c) { ... }
 * };
 */
template <typename T>
struct totally_ordered
    : less_than_comparable<T>, equality_comparable<T> {
};
