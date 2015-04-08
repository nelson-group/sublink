#pragma once
/** @file TreeUtil.hpp
 * @brief Useful functions for working with Subfind catalogs and
 *        SubLink merger trees.
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

#include <iostream>
#include <vector>

#include "../InputOutput/ReadTreeHDF5.hpp"
#include "TreeTypes.hpp"

/** Get some types from ReadTreeHDF5.hpp and make them our own. */
typedef Tree::Snapshot Snapshot;
typedef Tree::Subhalo Subhalo;

real_type periodic_dist(Subhalo sub1, Subhalo sub2, real_type box_size) {
  real_type d2 = 0;
  for (int i = 0; i < 3; ++i) {
    real_type dx = sub2.data().SubhaloPos[i] - sub1.data().SubhaloPos[i];
    if (std::abs(dx) > 0.5 * box_size)
      dx = dx - std::copysign(box_size, dx);
    d2 += dx*dx;
  }
  return std::sqrt(d2);
}

/** @brief Get the progenitor along the main branch at a given snapshot.
 * @param[in] sub The Subhalo of interest.
 * @param[in] snapnum The snapshot number of the progenitor of interest.
 * @pre @a sub is valid.
 * @return The latest progenitor along the main branch of @a sub such that
 *         @a result.snapnum() <= @a snapnum, or an invalid Subhalo if
 *         the main branch is truncated before reaching @a snapnum.
 */
Subhalo back_in_time(Subhalo sub, snapnum_type snapnum) {
  assert(sub.is_valid());
  while ((sub.snapnum() - snapnum > 0) && (sub.is_valid()))
    sub = sub.first_progenitor();
  return sub;
}

/** @brief Get the progenitors of @a primary and @a secondary at
 *         the moment of infall or earlier.
 * @param[in] primary The primary subhalo.
 * @param[in] secondary The secondary subhalo.
 * @return A pair with the main branch progenitors of @a primary
 *         and @a secondary at the latest snapshot such that
 *         @a primary and @a secondary do not belong to the same
 *         FoF group. Return a pair with invalid subhalos if one
 *         or both branches are truncated.
 *
 * @pre @a primary and @a secondary are valid subhalos.
 * @post new(@a primary).snapnum() == new(@a secondary).snapnum()
 */
std::pair<Subhalo, Subhalo> get_infall_pair(Subhalo primary,
    Subhalo secondary) {
  assert(primary.is_valid() && secondary.is_valid());

  while (true) {
    // If secondary branch is truncated, return invalid Subhalos.
    if (!secondary.is_valid())
      return std::make_pair(Subhalo(), Subhalo());

    // If primary branch is truncated, return invalid Subhalos.
    auto cur_snapnum = secondary.snapnum();
    primary = back_in_time(primary, cur_snapnum);
    if (!primary.is_valid())
      return std::make_pair(Subhalo(), Subhalo());

    // If primary skipped this snapshot, or is for some other reason found
    // at an earlier snapshot, "increment" secondary and try again.
    if (primary.snapnum() != cur_snapnum) {
      secondary = secondary.first_progenitor();
      continue;
    }
    // If we got here, both primary and secondary are valid Subhalos at
    // the same snapshot. If they do not belong to the same FoF group,
    // we have reached infall.
    if (  primary.first_subhalo_in_fof_group() !=
        secondary.first_subhalo_in_fof_group())
      return std::make_pair(primary, secondary);

    // Increment secondary
    secondary = secondary.first_progenitor();
  }
  assert(false);
}

/** @brief Get the progenitors of @a primary and @a secondary at
 *         the moment of the last virial crossing.
 * @param[in] primary The primary subhalo.
 * @param[in] secondary The secondary subhalo.
 * @return A pair with the main branch progenitors of @a primary
 *         and @a secondary at the latest snapshot such that
 *         @a primary and @a secondary are separated by more than
 *         @a R_{200,crit} of the parent halo of @a primary.
 *         Return a pair with invalid subhalos if one
 *         or both branches are truncated.
 *
 * @pre @a primary and @a secondary are valid subhalos.
 * @post new(@a primary).snapnum() == new(@a secondary).snapnum()
 */
std::pair<Subhalo, Subhalo> get_pair_last_crossing(Subhalo primary,
    Subhalo secondary, real_type box_size) {
  assert(primary.is_valid() && secondary.is_valid());

  while (true) {
    // If secondary branch is truncated, return invalid Subhalos.
    if (!secondary.is_valid())
      return std::make_pair(Subhalo(), Subhalo());

    // If primary branch is truncated, return invalid Subhalos.
    auto cur_snapnum = secondary.snapnum();
    primary = back_in_time(primary, cur_snapnum);
    if (!primary.is_valid())
      return std::make_pair(Subhalo(), Subhalo());

    // If primary skipped this snapshot, or is for some other reason found
    // at an earlier snapshot, "increment" secondary and try again.
    if (primary.snapnum() != cur_snapnum) {
      secondary = secondary.first_progenitor();
      continue;
    }
    // If we got here, both primary and secondary are valid Subhalos at
    // the same snapshot. If they are separated by more than R_{200,crit},
    // we have reached the moment of the last virial crossing.
    if (periodic_dist(primary, secondary, box_size) > primary.data().Group_R_Crit200)
      return std::make_pair(primary, secondary);

    // Increment secondary
    secondary = secondary.first_progenitor();
  }
  assert(false);
}

/** @brief Get the progenitors of @a primary and @a secondary at
 *         the moment of the first virial crossing.
 * @param[in] primary The primary subhalo.
 * @param[in] secondary The secondary subhalo.
 * @return A pair with the main branch progenitors of @a primary
 *         and @a secondary at the earliest snapshot such that
 *         @a primary and @a secondary are separated by more than
 *         @a R_{200,crit} of the parent halo of @a primary.
 *         Return a pair with invalid subhalos if one
 *         or both branches are truncated.
 *
 * @pre @a primary and @a secondary are valid subhalos.
 * @post new(@a primary).snapnum() == new(@a secondary).snapnum()
 */
std::pair<Subhalo, Subhalo> get_pair_first_crossing(Subhalo primary,
    Subhalo secondary, real_type box_size) {
  assert(primary.is_valid() && secondary.is_valid());

  // By default return invalid Subhalos.
  auto sub_pair = std::make_pair(Subhalo(), Subhalo());

  // Keep track of subhalo separation with respect to R200Crit.
  bool inside_r200 = true;

  while (true) {
    // If secondary branch is truncated, return invalid Subhalos.
    if (!secondary.is_valid())
      return std::make_pair(Subhalo(), Subhalo());

    // If primary branch is truncated, return invalid Subhalos.
    auto cur_snapnum = secondary.snapnum();
    primary = back_in_time(primary, cur_snapnum);
    if (!primary.is_valid())
      return std::make_pair(Subhalo(), Subhalo());

    // If primary skipped this snapshot, or is for some other reason found
    // at an earlier snapshot, "increment" secondary and try again.
    if (primary.snapnum() != cur_snapnum) {
      secondary = secondary.first_progenitor();
      continue;
    }

    // If we got here, both primary and secondary are valid Subhalos at
    // the same snapshot.
    if (periodic_dist(primary, secondary, box_size) > primary.data().Group_R_Crit200) {
      if (inside_r200) {
        // Overwrite info from later (time-wise) virial crossing
        sub_pair = std::make_pair(primary, secondary);
        inside_r200 = false;
      }
    }
    else {
      inside_r200 = true;
    }

    // Increment secondary
    secondary = secondary.first_progenitor();
  }

  return sub_pair;
}

/** @brief Get the progenitors of @a primary and @a secondary at
 *         the same snapshot as @a secondary or earlier.
 * @param[in] primary The primary subhalo.
 * @param[in] secondary The secondary subhalo.
 * @return A pair with the main branch progenitors of @a primary
 *         and @a secondary at the latest snapshot @a snapnum
 *         such that @a snapnum <= @a secondary.snapnum().
 *         Return a pair with invalid subhalos if one or both
 *         branches are truncated.
 *
 * @pre @a primary and @a secondary are valid subhalos.
 * @pre @a secondary.snapnum() <= @a primary.snapnum()
 * @post new(@a primary).snapnum() == new(@a secondary).snapnum()
 *
 * @note Implementation of this function is very similar to get_infall_pair().
 *       Both functions could rely on a more general function in which the two
 *       subhalos are "moved" back to a moment when a certain condition
 *       becomes true.
 */
std::pair<Subhalo, Subhalo> synchronize_subhalos(Subhalo primary,
    Subhalo secondary) {
  assert(primary.is_valid() && secondary.is_valid());
  assert(secondary.snapnum() <= primary.snapnum());

  while (true) {
    // If secondary branch is truncated, return invalid Subhalos.
    if (!secondary.is_valid())
      return std::make_pair(Subhalo(), Subhalo());

    // If primary branch is truncated, return invalid Subhalos.
    auto cur_snapnum = secondary.snapnum();
    primary = back_in_time(primary, cur_snapnum);
    if (!primary.is_valid())
      return std::make_pair(Subhalo(), Subhalo());

    // If primary skipped this snapshot, "increment" secondary and try again.
    if (primary.snapnum() != cur_snapnum) {
      secondary = secondary.first_progenitor();
      continue;
    }
    // If we got here, both primary and secondary are valid Subhalos at
    // the same snapshot.
    return std::make_pair(primary, secondary);
  }
  assert(false);
}

/** @brief Return true if @a secondary has already "infalled" into
 *         the same FoF group as @a primary; false otherwise.
 * @pre @a primary and @a secondary are valid subhalos.
 * @pre @a secondary.snapnum() <= @a primary.snapnum()
 */
bool after_infall(Subhalo primary, Subhalo secondary) {
  auto synced_pair = synchronize_subhalos(primary, secondary);
  primary = synced_pair.first;
  secondary = synced_pair.second;
  if (!primary.is_valid() || !secondary.is_valid())
    return false;
  assert(primary.first_subhalo_in_fof_group().is_valid());
  return (primary.first_subhalo_in_fof_group() ==
      secondary.first_subhalo_in_fof_group());
}

/** @brief Return true if @a prog lies along the main branch of @a desc.
 * @pre @a desc and @a prog are valid subhalos.
 */
bool along_main_branch(Subhalo desc, Subhalo prog) {
  return (prog.data().SubhaloID >= desc.data().SubhaloID) &&
         (prog.data().SubhaloID <= desc.data().MainLeafProgenitorID);
}

/** @brief Check if @a possible_desc is a descendant of @a possible_prog.
 * @pre @a desc and @a prog are valid subhalos.
 *
 * @note In this implementation a subhalo can be its own descendant.
 *       Maybe this function should be called "belongs_to_subtree".
 */
bool is_descendant(Subhalo desc, Subhalo prog) {
  return (prog.data().SubhaloID >= desc.data().SubhaloID) &&
         (prog.data().SubhaloID <= desc.data().LastProgenitorID);
}

/** Return the progenitor along the main branch which has the largest
 * mass (depending on merger tree type).
 * @pre @a sub is valid.
 */
Subhalo at_tmax(Subhalo sub) {
  auto sub_stmax = sub;
  sub = sub.first_progenitor();
  for ( ; sub.is_valid(); sub = sub.first_progenitor())
    if (sub.data().Mass > sub_stmax.data().Mass)
      sub_stmax = sub;
  return sub_stmax;
}

/** Return the progenitor along the main branch which has the largest
 * stellar mass.
 * @pre @a sub is valid.
 */
Subhalo at_stmax(Subhalo sub) {
  auto sub_stmax = sub;
  sub = sub.first_progenitor();
#ifdef COUNT_MERGERS
  for ( ; sub.is_valid(); sub = sub.first_progenitor())
    if (sub.data().SubhaloMassType[4] > sub_stmax.data().SubhaloMassType[4])
      sub_stmax = sub;
#else
  std::cerr << "Error: SubhaloMassType not included in trees." << std::endl;
  exit(1);
#endif

  return sub_stmax;
}

/** @brief Get the progenitors of @a primary and @a secondary at
 *         the moment when @a secondary reached its maximum stellar mass.
 * @param[in] primary The primary subhalo.
 * @param[in] secondary The secondary subhalo.
 * @return A pair with the main branch progenitors of @a primary
 *         and @a secondary at the moment when @a secondary reached
 *         its maximum stellar mass. The progenitor of @a primary
 *         may be found one snapshot earlier in case of skipping.
 *         Return a pair with invalid subhalos if the
 *         branch of @a primary was truncated.
 *
 * @pre @a primary and @a secondary are valid subhalos.
 */
std::pair<Subhalo, Subhalo> get_stmax_pair(Subhalo primary,
    Subhalo secondary) {
  assert(primary.is_valid() && secondary.is_valid());

  auto prog_stmax_2 = at_stmax(secondary);
  auto snapnum_stmax = prog_stmax_2.data().SnapNum;
  auto prog_stmax_1 = back_in_time(primary, snapnum_stmax);

  if (!prog_stmax_1.is_valid())
    return std::make_pair(Subhalo(), Subhalo());
  return std::make_pair(prog_stmax_1, prog_stmax_2);
}
