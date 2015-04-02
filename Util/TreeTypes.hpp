#pragma once
/** @file TreeTypes.hpp
 * @brief Define some common types used in many files.
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

/** @brief Type of particle IDs. */
typedef uint64_t part_id_type;
/** @brief Type of subhalo IDs in the merger trees. */
typedef int64_t sub_id_type;
/** @brief Type of merger tree IDs. */
typedef int64_t tree_id_type;
/** @brief Type of subhalo indices in the Subfind catalogs. */
typedef int32_t index_type;
/** @brief Type of subhalo lengths. */
typedef uint32_t sub_len_type;
/** @brief Type of snapshot numbers. */
typedef int16_t snapnum_type;
/** @brief Type of most physical quantities, e.g., masses. */
typedef float real_type;

/** @brief Type for representing arrays of floats (e.g., positions). */
template <int N = 1>
struct FloatArray {
  float elem[N];

  // Typedefs
  typedef float          value_type;
  typedef float&         reference;
  typedef const float&   const_reference;
  typedef float*         iterator;
  typedef const float*   const_iterator;

  // Accesors
  reference       operator[](std::size_t i)       { return elem[i]; }
  const_reference operator[](std::size_t i) const { return elem[i]; }

  /** Return the number of elements in this FloatArray. */
  static constexpr int size() {
    return N;
  }

};
