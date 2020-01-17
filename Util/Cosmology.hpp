#pragma once
/** @file Cosmology.hpp
 * @brief Some basic cosmology equations.
 *
 * @author Vicente Rodriguez-Gomez (v.rodriguez@irya.unam.mx)
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

/** @brief This class contains the cosmological parameters,
 * usually corresponding to a suite of cosmological simulations.
 */
class CosmologicalParameters {
public:
  // Main parameters:
  real_type Omega_m;
  real_type Omega_b;
  real_type Omega_L;
  real_type sigma_8;
  real_type h;
  real_type n_s;

  // Derived quantities:
  real_type a_mL;  // scale factor at radiation-matter equality
  real_type H0;  // Hubble constant in s^{-1}
  real_type H0_Gyr;  // Hubble constant in Gyr^-1
  const real_type H0_kpc_h = 0.1;  // Hubble constant in (km/s)/(kpc/h)

  /** @brief Constructor. Define cosmological parameters and calculate some
   * derived quantities.
   */
  CosmologicalParameters(const std::string& suite) {
    if (suite == "Illustris") {
      // WMAP-7, Komatsu et al. 2011 (Table 1, v2)
      Omega_m = 0.2726;
      Omega_b = 0.0456;
      Omega_L = 0.7274;
      sigma_8 = 0.809;
      h = 0.704;
      n_s = 0.963;
    }
    else if (suite == "IllustrisTNG") {
      // Planck 2015 XIII (Table 4, last column)
      Omega_m = 0.3089;
      Omega_b = 0.0486;
      Omega_L = 0.6911;
      sigma_8 = 0.8159;
      h = 0.6774;
      n_s = 0.9667;
    }
    else {
      std::cerr << "Unknown cosmological simulations: " << suite << std::endl;
      exit(1);
    }

    a_mL = std::pow(Omega_m / Omega_L, 1.0/3.0);
    H0 = h*100000.0/3.086e22;  // in s^{-1}
    H0_Gyr = H0 * (1e9*365.25*86400);  // in Gyr^-1
  }
};

/** @brief Return the Hubble parameter H(z) in units of s^{-1}. */
real_type H(const real_type z, const CosmologicalParameters& par) {
  return par.H0 * std::sqrt(par.Omega_m * std::pow(1.0+z, 3.0) + par.Omega_L);
}

/** @brief Return the Hubble parameter H(z) in units of Gyr^{-1}. */
real_type H_Gyr(const real_type z, const CosmologicalParameters& par) {
  return par.H0_Gyr * std::sqrt(par.Omega_m * std::pow(1.0+z, 3.0) + par.Omega_L);
}

/** @brief Return the Hubble parameter H(z) in units of (km/s)/(kpc/h). */
real_type H_kpc_h(const real_type z, const CosmologicalParameters& par) {
  return par.H0_kpc_h * std::sqrt(par.Omega_m * std::pow(1.0+z, 3.0) + par.Omega_L);
}

/** @brief Return the cosmic time (in s) at a given redshift (MvW, eq. 3.99).
 */
real_type t(const real_type z, const CosmologicalParameters& par) {
  return 1.0/par.H0 * 2.0/(3.0*sqrt(par.Omega_L)) * log((sqrt(par.Omega_L*pow(1.0+z, -3.0)) +
         sqrt(par.Omega_L*pow(1.0+z, -3.0) + par.Omega_m))/sqrt(par.Omega_m));
}

/** @brief Return cosmic time (in Gyr) at a given redshift.
 */
real_type t_Gyr(const real_type z, const CosmologicalParameters& par) {
  return t(z, par) / (1e9*365.25*86400);
}

/** @brief Read file with snapshot redshifts, one per line.
 * @return Vector with redshifts.
 *
 * We assume that the input file is located inside the working directory.
 */
std::vector<real_type> get_redshifts(const std::string& suite) {

  // Open file with redshifts
  std::stringstream ss;
  ss << "Redshifts" << suite << ".txt";
  std::string redshifts_filename = ss.str();
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
std::vector<real_type> get_times_Gyr(const std::string& suite) {
  auto redshifts_all = get_redshifts(suite);
  auto params = CosmologicalParameters(suite);
  std::vector<real_type> times_all;
  for (auto z : redshifts_all) {
    times_all.push_back(t_Gyr(z, params));
  }
  return times_all;
}

}  // end namespace cosmo

