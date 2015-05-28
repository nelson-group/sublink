#pragma once
/** @file SnapshotUtil.hpp
 * @brief Some functions for working with Arepo snapshot files and
 * Subfind catalogs.
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

#include <vector>
#include <fstream>
#include <cassert>
#include <cmath>
#include <algorithm>  // std::find

#include "../InputOutput/ReadSubfindHDF5.hpp"
#include "TreeTypes.hpp"

// Some constants
static constexpr real_type h = 0.704;
static constexpr real_type Omega_L = 0.7274;
static constexpr real_type Omega_m = 0.2726;
static constexpr real_type Omega_b = 0.0456;

// Some conversions
static constexpr real_type H0 = h*100000.0/3.086e22;  // in s^-1
static constexpr real_type H0_Gyr = H0 * (1e9*365.25*86400);  // in Gyr^-1

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

/** Create list of valid snapshots. */
std::vector<snapnum_type> get_valid_snapnums(
    const std::string& skipsnaps_filename,
    const snapnum_type& snapnum_first,
    const snapnum_type& snapnum_last) {

  // Open file with snapshot numbers, one per line
  std::ifstream infile (skipsnaps_filename.data());
  if (!infile.is_open()) {
    std::cerr << "Cannot open file: " << skipsnaps_filename << "\n";
    exit(1);
  }

  // Read snapshot numbers to be skipped
  std::vector<snapnum_type> invalid_snapnums;
  std::string line;
  while (infile >> line)
    invalid_snapnums.push_back(atoi(line.data()));
  infile.close();

  // Add valid snapshots to list
  std::vector<snapnum_type> valid_snapnums;
  for (auto snapnum = snapnum_first; snapnum <= snapnum_last; ++snapnum)
    if (std::find(invalid_snapnums.begin(), invalid_snapnums.end(), snapnum) ==
        invalid_snapnums.end())
      valid_snapnums.push_back(snapnum);

  return valid_snapnums;
}

/** @brief Read file with snapshot redshifts, one per line.
 * @return Vector with redshifts.
 *
 * We assume that the input file is located inside the working directory.
 */
std::vector<real_type> get_redshifts() {

  // Open file with redshifts
  std::string redshifts_filename = "redshifts_illustris.txt";
  std::ifstream infile(redshifts_filename.data());
  if (!infile.is_open()) {
    std::cerr << "Unable to open file: " << redshifts_filename << std::endl;
    exit(1);
  }

  // Read redshifts
  std::vector<real_type> redshifts_all;
  std::string line;
  while (infile >> line)
    redshifts_all.push_back(atof(line.data()));
  infile.close();

  return redshifts_all;
}

/** @brief Return the cosmic time for a given redshift.
 *
 * Equation from MvW (Eq. 3.99)
 */
real_type t(const real_type z) {
  return 1.0/H0 * 2.0/(3.0*sqrt(Omega_L)) * log((sqrt(Omega_L*pow(1.0+z, -3.0)) +
         sqrt(Omega_L*pow(1.0+z, -3.0) + Omega_m))/sqrt(Omega_m));
}

/** Return cosmic time for a given redshift in Gyr. */
real_type t_Gyr(const real_type z) {
  return t(z) / (1e9*365.25*86400);
}

/** Get cosmic time in Gyr for each snapshot. */
std::vector<real_type> get_times_Gyr() {
  auto redshifts_all = get_redshifts();
  std::vector<real_type> times_all;
  for (auto z : redshifts_all) {
    times_all.push_back(t_Gyr(z));
  }
  return times_all;
}
