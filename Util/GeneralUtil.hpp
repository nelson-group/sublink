#pragma once
/** @file GeneralUtil.hpp
 * @brief General utilities, some of them taken from CS207.
 */

#include <iostream>
#include <vector>
#include <chrono>   // Wall clock time
#include <ctime>    // CPU time

/** @brief Class for measuring wall clock time. */
class WallClock {
public:
  /** Constructor. */
  WallClock() : start_(std::chrono::system_clock::now()) {
  }
  /** Start this clock. */
  void start() {
    start_ = std::chrono::system_clock::now();
  }
  /** Return elapsed time in seconds. */
  double seconds() {
    return std::chrono::duration<double>(
        std::chrono::system_clock::now()-start_).count();
  }
private:
  std::chrono::system_clock::time_point start_;
};

/** @brief Class for measuring CPU time. */
class CPUClock {
public:
  /** Constructor. */
  CPUClock() : start_(std::clock()) {
  }
  /** Start this clock. */
  void start() {
    start_ = std::clock();
  }
  /** Return elapsed time in seconds. */
  double seconds() {
    return 1.0 * (std::clock()-start_) / CLOCKS_PER_SEC;
  }
private:
  std::clock_t start_;
};

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
