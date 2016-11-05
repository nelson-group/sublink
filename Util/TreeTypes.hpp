#pragma once
/** @file TreeTypes.hpp
 * @brief Define some common types used in many files.
 *
 * @author Vicente Rodriguez-Gomez (vrg@jhu.edu)
 */

#include <array> // (C++11)

/** @brief Type of particle IDs. MUST BE uint32_t for non-LONGIDS runs. */
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

// Replace old FloatArray and DoubleArray types with std::array
// using alias declarations (C++11):

/** @brief Type for representing arrays of floats (e.g., positions). */
template <int N>
using FloatArray = std::array<float, N>;

/** @brief Type for representing arrays of doubles (e.g., positions). */
template <int N>
using DoubleArray = std::array<double, N>;
