#pragma once
/** @file ParticleMatcher.hpp
 * @brief Define a class for matching particles across different snapshots.
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>  // setfill, setw, setprecision
#include <chrono>   // Wall clock time
#include <ctime>    // CPU time
#include <cmath>    // pow
#include <cassert>

#include "../InputOutput/ReadArepoHDF5.hpp"
#include "../InputOutput/ReadSubfindHDF5.hpp"
#include "../Util/SnapshotUtil.hpp"
#include "../Util/GeneralUtil.hpp"

// Determines how important is the contribution from the innermost
// particles in a subhalo when finding a descendant.
static constexpr float alpha_weight = -1;

class ParticleMatcher {

public:
  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  /** @brief Type of this ParticleMatcher. */
  typedef ParticleMatcher particle_matcher_type;

  /** @brief Type of particle IDs. */
  typedef uint64_t part_id_type;
  /** @brief Type of subhalo indices in the Subfind catalogs. */
  typedef int32_t index_type;
  /** @brief Type of subhalo lengths. */
  typedef uint32_t sub_len_type;
  /** @brief Type of snapshot numbers. */
  typedef int16_t snapnum_type;
  /** @brief Type of most physical quantities, e.g., masses. */
  typedef float real_type;

  /** Predeclare Snapshot type. */
  class Snapshot;
  /** @brief Synonym for Snapshot. */
  typedef Snapshot snapshot_type;

  /** Datatype for simulation particles. */
  struct ParticleInfo{
    part_id_type id;
    index_type sub_index;
    real_type weight;
    /** Constructor. */
    ParticleInfo(part_id_type id_, index_type sub_index_, real_type weight_)
        : id(id_), sub_index(sub_index_), weight(weight_) {
    }
  };

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  /** Constructor. */
  ParticleMatcher(const std::string& basedir1, const std::string& basedir2,
      const snapnum_type snapnum1, const snapnum_type snapnum2,
      const std::string& tracking_scheme)
      : snap1_(nullptr), snap2_(nullptr), data_() {

    // Create Snapshot objects.
    snap1_ = new Snapshot(this, basedir1, snapnum1, tracking_scheme);
    if (snapnum2 != -1) {
      snap2_ = new Snapshot(this, basedir2, snapnum2, tracking_scheme);
    }

  }
  /** Destructor. */
  ~ParticleMatcher() {
    delete snap1_;
    if (snap2_ != nullptr)
      delete snap2_;
  }


  ///////////////
  // SNAPSHOTS //
  ///////////////

  class Snapshot {


  public:


  private:
    friend class ParticleMatcher;
    // Pointer back to the ParticleMatcher container.
    ParticleMatcher* pm_;
    // Number of particles in each subhalo.
    std::vector<uint32_t> sub_len_;
    // Mass of each subhalo.
    std::vector<real_type> sub_mass_;
    // Index of descendant subhalo.
    std::vector<index_type> descendants_;
    // Score of unique descendant.
    std::vector<real_type> max_scores_;

    /** Private constructor. */
    Snapshot(const ParticleMatcher* pm, const std::string& basedir,
        const snapnum_type snapnum, const std::string& tracking_scheme)
        : pm_(const_cast<ParticleMatcher*>(pm)),
          sub_len_(),
          sub_mass_(),
          descendants_(),
          max_scores_() {

      // Read data
      read_ids(basedir, snapnum, tracking_scheme);

    }


    void read_ids(const std::string& basedir, const snapnum_type snapnum,
        const std::string& tracking_scheme) {
      // For performance checks
      WallClock clock_all;
      WallClock clock;

      // Define particle types
      std::vector<int> parttypes;
      if (tracking_scheme == "Subhalos")
        parttypes = {1};  // DM
      else if (tracking_scheme == "Galaxies")
        parttypes = {0, 4};  // gas and stars
      else
        assert(false);
      unsigned num_parttypes = parttypes.size();

      // Load some subhalo info and initialize member variables
      std::cout << "Loading subhalo info...\n";
      clock.start();
      auto nsubs = subfind::get_scalar_attribute<uint32_t>(
          basedir, snapnum, "Nsubgroups_Total");
      std::vector<std::vector<uint32_t>> sub_len_parttype;
      std::vector<std::vector<uint32_t>> sub_offset_parttype;
      for (unsigned l = 0; l < num_parttypes; ++l) {
        sub_len_parttype.push_back(subfind::read_block<uint32_t>(
            basedir, snapnum, "Subhalo", "SubhaloLenType", parttypes[l]));
        sub_offset_parttype.push_back(calculate_subhalo_offsets(
            basedir, snapnum, parttypes[l]));
      }
      // Initialize member variables
      if (tracking_scheme == "Subhalos") {
        sub_len_ = sub_len_parttype[0];  // DM
        sub_mass_ = subfind::read_block<real_type>(
            basedir, snapnum, "Subhalo", "SubhaloMassType", parttypes[0]);
      }
      else {  // Galaxies
        sub_len_ = std::vector<uint32_t>(nsubs, 0);
        sub_mass_ = std::vector<real_type>(nsubs, 0);
      }
      descendants_ = std::vector<index_type>(nsubs, -1);
      max_scores_ = std::vector<real_type>(nsubs, 0);
      std::cout << "Time: " << clock.seconds() << " s.\n";

      // Number of particles in subhalos (i.e., exclude unbound particles)
      uint64_t npart = 0;
      for (unsigned l = 0; l < num_parttypes; ++l)
        for (uint32_t i = 0; i < nsubs; ++i)
          npart += sub_len_parttype[l][i];

      // Reserve space in memory
      uint64_t data_start = pm_->data_.size();
      pm_->data_.reserve(data_start + npart);

      for (unsigned l = 0; l < num_parttypes; ++l) {
        // Load particle IDs
        std::cout << "Loading particle IDs...\n";
        clock.start();
        auto part_id = arepo::read_block<part_id_type>(
                  basedir, snapnum, "ParticleIDs", parttypes[l]);
        std::cout << "Time: " << clock.seconds() << " s.\n";

        // Associate particles with subhalos
        std::cout << "Associating particles with subhalos...\n";
        clock.start();
        if (parttypes[l] == 1) {  // DM
          for (uint32_t sub_uindex = 0; sub_uindex < nsubs; ++sub_uindex) {
            part_id_type snap_count = sub_offset_parttype[l][sub_uindex];
            for (uint32_t i = 0; i < sub_len_parttype[l][sub_uindex]; ++i) {
              pm_->data_.emplace_back(ParticleInfo(
                  part_id[snap_count],
                  sub_uindex,
                  std::pow(static_cast<real_type>(i+1), alpha_weight)));
              ++snap_count;
            }
          }
        }
        else {  // stars or gas
          auto part_mass = arepo::read_block<real_type>(
                      basedir, snapnum, "Masses", parttypes[l]);
          for (uint32_t sub_uindex = 0; sub_uindex < nsubs; ++sub_uindex) {
            part_id_type snap_count = sub_offset_parttype[l][sub_uindex];
            for (uint32_t i = 0; i < sub_len_parttype[l][sub_uindex]; ++i) {
              pm_->data_.emplace_back(ParticleInfo(
                  part_id[snap_count],
                  sub_uindex,
                  part_mass[snap_count] * std::pow(
                      static_cast<real_type>(i+1), alpha_weight)));
              sub_len_[sub_uindex] += 1;
              sub_mass_[sub_uindex] += part_mass[snap_count];
              ++snap_count;
            }
          }
        }
        std::cout << "Time: " << clock.seconds() << " s.\n";
        std::cout << "Finished for parttype " << parttypes[l] << ".\n";
      }
      assert(data_start + npart == pm_->data_.size());
      std::cout << "Finished reading snapshot " << snapnum << ".\n";
      std::cout << "Total time: " << clock_all.seconds() << " s.\n\n";
    }
  };


private:

  //////////////////////////////
  // PRIVATE MEMBER FUNCTIONS //
  //////////////////////////////



  //////////////////////////////
  // PRIVATE MEMBER VARIABLES //
  //////////////////////////////

  Snapshot* snap1_;
  Snapshot* snap2_;
  std::vector<ParticleInfo> data_;

};
