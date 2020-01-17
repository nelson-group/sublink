#pragma once
/** @file HybridSort.hpp
 * @brief Implement a new sorting function with performance somewhere in between
 *        the serial and parallel versions of std::stable_sort.
 *
 * @author Vicente Rodriguez-Gomez (v.rodriguez@irya.unam.mx)
 */

#include <iostream>
#include <cassert>
#include <algorithm> // merge
#include <parallel/algorithm> // parallel stable_sort

#include "TreeTypes.hpp"
#include "GeneralUtil.hpp"

/** @brief A merge sort-like algorithm with performance somewhere in between
 *         the serial and parallel versions of std::stable_sort.
 */
template <typename RandomAccessIterator, typename Compare>
void hybrid_sort(RandomAccessIterator first, RandomAccessIterator last,
    Compare comp) {
  assert(last - first > 1);
  RandomAccessIterator middle = first + (last - first) / 2;
  __gnu_parallel::stable_sort(first, middle, comp);
  __gnu_parallel::stable_sort(middle, last, comp);
  std::inplace_merge(first, middle, last, comp);
}

template <typename RandomAccessIterator>
void hybrid_sort(RandomAccessIterator first, RandomAccessIterator last) {
  hybrid_sort(first, last,
      std::less<typename std::iterator_traits<RandomAccessIterator>::value_type>());
}
