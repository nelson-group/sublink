#pragma once
/** @file SnapshotUtil.hpp
 * @brief Some functions for working with Arepo snapshot files and
 * Subfind catalogs.
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

#include <vector>
#include <cassert>

#include "TreeTypes.hpp"

/** @brief Function to calculate subhalo offsets.
 *
 * In this context, the subhalo offset is the index of the
 * first particle of a given type that belongs to each subhalo.
 */
std::vector<uint32_t> calculate_subhalo_offsets(const std::string& basedir,
    const snapnum_type snapnum, const int parttype) {

  // Load some FoF group and subhalo info
  auto group_nsubs = subfind::read_block<uint32_t>(basedir, snapnum, "Group",
      "GroupNsubs", -1);
  auto ngroups = subfind::get_scalar_attribute<uint32_t>(basedir, snapnum,
      "Ngroups_Total");
  auto nsubs = subfind::get_scalar_attribute<uint32_t>(basedir, snapnum,
      "Nsubgroups_Total");
  auto group_len = subfind::read_block<uint32_t>(basedir, snapnum, "Group",
      "GroupLenType", parttype);
  auto sub_len = subfind::read_block<uint32_t>(basedir, snapnum, "Subhalo",
      "SubhaloLenType", parttype);
  std::vector<uint32_t> group_offset(ngroups, 0);
  std::vector<uint32_t> sub_offset(nsubs, 0);

  // Calculate offsets
  uint32_t k = 0;
  for (uint32_t i = 0; i < ngroups; ++i) {
    if (i>0)
      group_offset[i] = group_offset[i-1] + group_len[i-1];
    if (group_nsubs[i] > 0) {
      sub_offset[k] = group_offset[i];
      ++k;
      for (uint32_t j = 1; j < group_nsubs[i]; ++j) {
        sub_offset[k] = sub_offset[k-1] + sub_len[k-1];
        ++k;
      }
    }
  }
  assert(k == nsubs);

  return sub_offset;
}
