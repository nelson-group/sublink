#pragma once
/** @file Cosmology.hpp
 * @brief Some basic cosmology equations.
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

#include <vector>
#include <fstream>
#include <cassert>
#include <cmath>

#include "TreeTypes.hpp"

/** @namespace cosmo
 * @brief Namespace for basic cosmological constants and functions.
 */
namespace cosmo {

// Some constants
static constexpr real_type h = 0.704;
static constexpr real_type Omega_L = 0.7274;
static constexpr real_type Omega_m = 0.2726;
static constexpr real_type Omega_b = 0.0456;

// Some conversions
static constexpr real_type H0 = h*100000.0/3.086e22;  // in s^-1
static constexpr real_type H0_Gyr = H0 * (1e9*365.25*86400);  // in Gyr^-1
static constexpr real_type H0_kpc_h = 0.1;    // in (km/s)/(kpc/h)

/** @brief Return the Hubble parameter H(z) in units of s^{-1}. */
real_type H(const real_type z) {
  return H0 * std::sqrt(Omega_m * std::pow(1.0+z, 3.0) + Omega_L);
}

/** @brief Return the Hubble parameter H(z) in units of Gyr^{-1}. */
real_type H_Gyr(const real_type z) {
  return H0_Gyr * std::sqrt(Omega_m * std::pow(1.0+z, 3.0) + Omega_L);
}

/** @brief Return the Hubble parameter H(z) in units of (km/s)/(kpc/h). */
real_type H_kpc_h(const real_type z) {
  return H0_kpc_h * std::sqrt(Omega_m * std::pow(1.0+z, 3.0) + Omega_L);
}

/** @brief Return the cosmic time for a given redshift.
 *
 * Equation from MvW (Eq. 3.99).
 */
real_type t(const real_type z) {
  return 1.0/H0 * 2.0/(3.0*sqrt(Omega_L)) * log((sqrt(Omega_L*pow(1.0+z, -3.0)) +
         sqrt(Omega_L*pow(1.0+z, -3.0) + Omega_m))/sqrt(Omega_m));
}

/** @brief Return cosmic time for a given redshift in Gyr.
 */
real_type t_Gyr(const real_type z) {
  return t(z) / (1e9*365.25*86400);
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

/** Get cosmic time in Gyr for each snapshot. */
std::vector<real_type> get_times_Gyr() {
  auto redshifts_all = get_redshifts();
  std::vector<real_type> times_all;
  for (auto z : redshifts_all) {
    times_all.push_back(t_Gyr(z));
  }
  return times_all;
}

}  // end namespace cosmo

