#pragma once
/** @file ParticleMatcher.hpp
 * @brief Define a class for matching particles across different snapshots.
 *
 * @author Vicente Rodriguez-Gomez (vrg@jhu.edu)
 */
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>  // setfill, setw, setprecision
#include <cmath>    // pow
#include <cassert>

#include <algorithm>  // find, stable_sort
#ifdef USE_OPENMP
#include <parallel/algorithm>  // parallel stable_sort
#include "../Util/HybridSort.hpp"  // "hybrid" sort
#endif

#include "../InputOutput/ReadArepoHDF5.hpp"
#include "../InputOutput/ReadSubfindHDF5.hpp"
#include "../InputOutput/GeneralHDF5.hpp"
#include "../Util/SnapshotUtil.hpp"
#include "../Util/GeneralUtil.hpp"
#include "../Util/TreeTypes.hpp"

// Determines how important is the contribution from the innermost
// particles in a subhalo when finding a descendant.
static constexpr float alpha_weight = -1;

//////////////////////
// TYPE DEFINITIONS //
//////////////////////

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

/** Comparison function to sort by particle ID. */
bool compareByID(const ParticleInfo& a, const ParticleInfo& b) {
  return a.id < b.id;
}

/** Structure to represent descendant candidates. */
struct Candidate {
  /** Constructor. */
  Candidate(index_type sub_index, real_type weight)
      : value_(std::make_pair(sub_index, weight)) {
  }

  /** Return the Subfind ID of this descendant candidate. */
  index_type index() const {
    return value_.first;
  }
  /** Return the score of this descendant candidate. */
  real_type score() const {
    return value_.second;
  }
  /** Add contribution from a particle. */
  void add_to_score(const real_type& weight) {
    value_.second += weight;
  }

  /** Custom equality operator for use with std::find */
  bool operator==(const index_type& other_index) const {
    return value_.first == other_index;
  }
  /** Custom less-than operator for use with std::max_element. */
  bool operator<(const Candidate& other) const {
    return value_.second < other.value_.second;
  }
private:
  // Representation: a <subfind_id, score> pair
  std::pair<index_type, real_type> value_;
};

////////////////////////////
// PARTICLE MATCHER CLASS //
////////////////////////////

/** @class ParticleMatcher
 * @brief Class for matching particles between two different snapshots.
 */
class ParticleMatcher {
public:
  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  /** @brief Type of this ParticleMatcher. */
  typedef ParticleMatcher particle_matcher_type;

  /** Predeclare Snapshot type. */
  class Snapshot;
  /** @brief Synonym for Snapshot. */
  typedef Snapshot snapshot_type;

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
      match_particles();
    }
  }

  /** Destructor. */
  ~ParticleMatcher() {
    delete snap1_;
    if (snap2_ != nullptr)
      delete snap2_;
  }

  /////////////////////////////
  // PUBLIC MEMBER FUNCTIONS //
  /////////////////////////////

  /** @brief Write to an HDF5 file. */
  void write_to_file(const std::string& writepath, bool writemisc = true) const {
    snap1_->write_to_file(writepath, writemisc);
  }

  ///////////////
  // SNAPSHOTS //
  ///////////////

  /** @class Snapshot
   * @brief Class to represent snapshots.
   */
  class Snapshot {
  public:
    /** Default constructor. Creates invalid Snapshot. */
    Snapshot() : pm_(nullptr), basedir_(), snapnum_(-1), sub_len_(),
        sub_mass_(), sub_grnr_(), descendants_(), first_scores_(),
        second_scores_() {
    }
    /** Default destructor. */
    ~Snapshot() = default;

    /** Return the number of subhalos in this Snapshot. */
    uint32_t nsubs() const {
      return sub_len_.size();
    }

  private:
    friend class ParticleMatcher;
    // Pointer back to the ParticleMatcher container.
    ParticleMatcher* pm_;
    // Directory containing the snapshot files.
    std::string basedir_;
    // Snapshot number.
    snapnum_type snapnum_;
    // Number of particles in each subhalo.
    std::vector<uint32_t> sub_len_;
    // Mass of each subhalo.
    std::vector<real_type> sub_mass_;
    // Index of parent FoF group
    std::vector<uint32_t> sub_grnr_;
    // Index of descendant subhalo.
    std::vector<index_type> descendants_;
    // Score of unique descendant.
    std::vector<real_type> first_scores_;
    // Score of second-best descendant candidate.
    std::vector<real_type> second_scores_;

    /** Private constructor. */
    Snapshot(const ParticleMatcher* pm, const std::string& basedir,
        const snapnum_type snapnum, const std::string& tracking_scheme)
        : pm_(const_cast<ParticleMatcher*>(pm)), basedir_(basedir),
          snapnum_(snapnum), sub_len_(), sub_mass_(), sub_grnr_(),
          descendants_(), first_scores_(), second_scores_() {
      // Read data
      read_ids(tracking_scheme);
    }

    /** @brief Read particle IDs and other information. */
    void read_ids(const std::string& tracking_scheme) {
      // For performance checks
      WallClock wall_clock_all;
      WallClock wall_clock;

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
      wall_clock.start();
      auto nsubs = subfind::get_scalar_attribute<uint32_t>(
          basedir_, snapnum_, "Nsubgroups_Total");
      std::vector<std::vector<uint32_t>> sub_len_parttype;
      std::vector<std::vector<uint64_t>> sub_offset_parttype;
      for (unsigned l = 0; l < num_parttypes; ++l) {
        sub_len_parttype.push_back(subfind::read_block<uint32_t>(
            basedir_, snapnum_, "Subhalo", "SubhaloLenType", parttypes[l]));
        sub_offset_parttype.push_back(calculate_subhalo_offsets(
            basedir_, snapnum_, parttypes[l]));
      }
      // Initialize member variables
      if (tracking_scheme == "Subhalos") {
        sub_len_ = sub_len_parttype[0];  // DM
        sub_mass_ = subfind::read_block<real_type>(
            basedir_, snapnum_, "Subhalo", "SubhaloMassType", parttypes[0]);
      }
      else {  // Galaxies
        sub_len_ = std::vector<uint32_t>(nsubs, 0);
        sub_mass_ = std::vector<real_type>(nsubs, 0);
      }
      sub_grnr_ = subfind::read_block<uint32_t>(
          basedir_, snapnum_, "Subhalo", "SubhaloGrNr", -1);
      descendants_ = std::vector<index_type>(nsubs, -1);
      first_scores_  = std::vector<real_type>(nsubs, 0);
      second_scores_ = std::vector<real_type>(nsubs, 0);
      std::cout << "Time: " << wall_clock.seconds() << " s.\n";

      // Reserve space in memory when dealing with DM
      if (tracking_scheme == "Subhalos") {
        part_id_type npart_in_subhalos = 0;
        for (uint32_t i = 0; i < nsubs; ++i)
          npart_in_subhalos += sub_len_parttype[0][i];
        pm_->data_.reserve(pm_->data_.size() + npart_in_subhalos);
      }

      for (unsigned l = 0; l < num_parttypes; ++l) {
        // Load particle IDs
        std::cout << "Loading particle IDs...\n";
        wall_clock.start();
        uint64_t nread = sub_offset_parttype[l][nsubs] + 1;

        auto part_id = arepo::read_block<part_id_type>(
                  basedir_, snapnum_, "ParticleIDs", parttypes[l], nread);
        std::cout << "Time: " << wall_clock.seconds() << " s.\n";

        // Associate particles with subhalos
        std::cout << "Associating particles with subhalos...\n";
        wall_clock.start();
        if (parttypes[l] == 1) {  // DM
          for (uint32_t sub_uindex = 0; sub_uindex < nsubs; ++sub_uindex) {
            part_id_type snap_count = sub_offset_parttype[l][sub_uindex];
            for (uint32_t i = 0; i < sub_len_parttype[l][sub_uindex]; ++i) {
              pm_->data_.emplace_back(
                  part_id[snap_count],
                  sub_uindex,
                  std::pow(static_cast<real_type>(i+1), alpha_weight));
              ++snap_count;
            }
          }
        }
        else if (parttypes[l] == 0) {  // gas
          auto part_mass = arepo::read_block<real_type>(
                      basedir_, snapnum_, "Masses", parttypes[l], nread);
          auto part_sfr = arepo::read_block<real_type>(
                      basedir_, snapnum_, "StarFormationRate", parttypes[l], nread);
          for (uint32_t sub_uindex = 0; sub_uindex < nsubs; ++sub_uindex) {
            part_id_type snap_count = sub_offset_parttype[l][sub_uindex];
            for (uint32_t i = 0; i < sub_len_parttype[l][sub_uindex]; ++i) {
              // Only consider star-forming elements
              if (part_sfr[snap_count] > 0) {
                pm_->data_.emplace_back(
                    part_id[snap_count],
                    sub_uindex,
                    part_mass[snap_count] * std::pow(
                        static_cast<real_type>(i+1), alpha_weight));
                sub_len_[sub_uindex] += 1;
                sub_mass_[sub_uindex] += part_mass[snap_count];
              }
              ++snap_count;
            }
          }
        }
        else if (parttypes[l] == 4) {  // stars
          auto part_mass = arepo::read_block<real_type>(
                      basedir_, snapnum_, "Masses", parttypes[l], nread);
          for (uint32_t sub_uindex = 0; sub_uindex < nsubs; ++sub_uindex) {
            part_id_type snap_count = sub_offset_parttype[l][sub_uindex];
            for (uint32_t i = 0; i < sub_len_parttype[l][sub_uindex]; ++i) {
              pm_->data_.emplace_back(
                  part_id[snap_count],
                  sub_uindex,
                  part_mass[snap_count] * std::pow(
                      static_cast<real_type>(i+1), alpha_weight));
              sub_len_[sub_uindex] += 1;
              sub_mass_[sub_uindex] += part_mass[snap_count];
              ++snap_count;
            }
          }
        }
        else {
          assert(false);
        }
        std::cout << "Time: " << wall_clock.seconds() << " s.\n";
        std::cout << "Finished for parttype " << parttypes[l] << ".\n";
      }
      std::cout << "Finished reading snapshot " << snapnum_ << ".\n";
      std::cout << "Total time: " << wall_clock_all.seconds() << " s.\n\n";
    }

    /** @brief Write to an HDF5 file. */
    void write_to_file(const std::string& writepath, bool writemisc) const {
      // Create filename
      std::stringstream tmp_stream;
      tmp_stream << writepath << "_" <<
          std::setfill('0') << std::setw(3) << snapnum_ << ".hdf5";
      std::string writefilename = tmp_stream.str();

      // Write to file
      std::cout << "Writing to file...\n";
      WallClock wall_clock;
      H5::H5File file(writefilename, H5F_ACC_TRUNC);
      if (writemisc) {
        add_array(file, sub_len_, "SubhaloLen", H5::PredType::NATIVE_UINT32);
        add_array(file, sub_mass_, "SubhaloMass", H5::PredType::NATIVE_FLOAT);
        add_array(file, sub_grnr_, "SubhaloGrNr", H5::PredType::NATIVE_UINT32);
      }
      add_array(file, descendants_, "DescendantIndex", H5::PredType::NATIVE_INT32);
      add_array(file, first_scores_, "FirstScore", H5::PredType::NATIVE_FLOAT);
      add_array(file, second_scores_, "SecondScore", H5::PredType::NATIVE_FLOAT);
      file.close();
      std::cout << "Time: " << wall_clock.seconds() << " s.\n";
    }
  };

private:

  //////////////////////////////
  // PRIVATE MEMBER VARIABLES //
  //////////////////////////////

  Snapshot* snap1_;
  Snapshot* snap2_;
  std::vector<ParticleInfo> data_;

  //////////////////////////////
  // PRIVATE MEMBER FUNCTIONS //
  //////////////////////////////

  /** @brief Match particles between the two snapshots. */
  void match_particles() {
    // Sort array
    std::cout << "Sorting array...\n";
    WallClock wall_clock;
    CPUClock cpu_clock;
#ifdef USE_OPENMP
//    __gnu_parallel::stable_sort(data_.begin(), data_.end(), compareByID);
    hybrid_sort(data_.begin(), data_.end(), compareByID);
#else
    std::stable_sort(data_.begin(), data_.end(), compareByID);
#endif
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";
    std::cout << "Speedup: " << cpu_clock.seconds()/wall_clock.seconds() << ".\n";

    // Iterate over particle coincidences to calculate scores
    std::cout << "Calculating scores...\n";
    wall_clock.start();
    std::vector<std::vector<Candidate>> scores(snap1_->nsubs());
    for (auto data_it = data_.begin(); data_it+1 < data_.end(); ++data_it) {

      // Check for bad part_id_type (see 2016/04/25 commit)
      assert(data_it->id != 0);

      // Only care about repeated IDs
      if (data_it->id != (data_it+1)->id)
        continue;

      // Add to score of current descendant candidate
      index_type sub_index1 = data_it->sub_index;  // progenitor
      index_type sub_index2 = (data_it+1)->sub_index;  // candidate

      auto& cur_cands = scores[sub_index1];
      auto it = std::find(cur_cands.begin(), cur_cands.end(), sub_index2);
      if (it == cur_cands.end())
        cur_cands.emplace_back(sub_index2, data_it->weight);
      else
        it->add_to_score(data_it->weight);

      // Sanity checks
      assert(sub_index1 < (int)snap1_->nsubs());
      assert(sub_index2 < (int)snap2_->nsubs());
      if ( (data_it+1)->id == (data_it+2)->id )
        std::cout << "WARNING DUPLICATE ID: it+1 id=" << (data_it+1)->id << " sub=" << (data_it+1)->sub_index 
                  << " it+2 id=" << (data_it+2)->id << " sub=" << (data_it+2)->sub_index << std::endl;
        //assert( (data_it+1)->id != (data_it+2)->id ); // not true for L75n1820TNG due to duplicate PT4 IDs
    }
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    // Determine descendants
    std::cout << "Determining descendants...\n";
    wall_clock.start();
    for (uint32_t sub_index1 = 0; sub_index1 < snap1_->nsubs(); ++sub_index1) {
      auto& cur_cands = scores[sub_index1];

      if (cur_cands.size() == 0)
        continue;

      real_type first_score = 0;
      real_type second_score = 0;
      index_type desc_index = -1;
      for (auto it = cur_cands.begin(); it != cur_cands.end(); ++it) {
        auto cur_score = it->score();
        if (cur_score > first_score) {
          second_score = first_score;
          first_score = cur_score;
          desc_index = it->index();
        }
        else if (cur_score > second_score) {
          second_score = cur_score;
        }
      }
      snap1_->descendants_[sub_index1] = desc_index;
      snap1_->first_scores_[sub_index1] = first_score;
      snap1_->second_scores_[sub_index1] = second_score;
    }
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  }
};
