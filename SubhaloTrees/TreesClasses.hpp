#pragma once
/** @file TreesClasses.hpp
 * @brief Define classes for constructing SubLink merger trees.
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cassert>

#include <algorithm>  // stable_sort, lower_bound
#ifdef USE_OPENMP
#include <parallel/algorithm>  // parallel stable_sort
#endif

#include "../InputOutput/GeneralHDF5.hpp"
#include "../Util/SnapshotUtil.hpp"
#include "../Util/GeneralUtil.hpp"
#include "../Util/TreeTypes.hpp"

class AllTrees {
public:
  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  /** @brief Format of the subhalo data, consisting of the "minimal" fields.
   */
  struct DataFormat {
    // Fields corresponding to "minimal" data format.
    sub_id_type SubhaloID;
    sub_id_type SubhaloIDRaw;
    sub_id_type LastProgenitorID;
    sub_id_type MainLeafProgenitorID;
    sub_id_type RootDescendantID;
    tree_id_type TreeID;
    snapnum_type SnapNum;
    sub_id_type FirstProgenitorID;
    sub_id_type NextProgenitorID;
    sub_id_type DescendantID;
    sub_id_type FirstSubhaloInFOFGroupID;
    sub_id_type NextSubhaloInFOFGroupID;
    sub_len_type NumParticles;
    real_type Mass;
    real_type MassHistory;
    index_type SubfindID;
    /** Constructor. */
    DataFormat(sub_id_type SubhaloID_,
                sub_id_type SubhaloIDRaw_,
                sub_id_type LastProgenitorID_,
                sub_id_type MainLeafProgenitorID_,
                sub_id_type RootDescendantID_,
                tree_id_type TreeID_,
                snapnum_type SnapNum_,
                sub_id_type FirstProgenitorID_,
                sub_id_type NextProgenitorID_,
                sub_id_type DescendantID_,
                sub_id_type FirstSubhaloInFOFGroupID_,
                sub_id_type NextSubhaloInFOFGroupID_,
                sub_len_type NumParticles_,
                real_type Mass_,
                real_type MassHistory_,
                index_type SubfindID_)
        : SubhaloID(SubhaloID_),
          SubhaloIDRaw(SubhaloIDRaw_),
          LastProgenitorID(LastProgenitorID_),
          MainLeafProgenitorID(MainLeafProgenitorID_),
          RootDescendantID(RootDescendantID_),
          TreeID(TreeID_),
          SnapNum(SnapNum_),
          FirstProgenitorID(FirstProgenitorID_),
          NextProgenitorID(NextProgenitorID_),
          DescendantID(DescendantID_),
          FirstSubhaloInFOFGroupID(FirstSubhaloInFOFGroupID_),
          NextSubhaloInFOFGroupID(NextSubhaloInFOFGroupID_),
          NumParticles(NumParticles_),
          Mass(Mass_),
          MassHistory(MassHistory_),
          SubfindID(SubfindID_) {
    }
  };

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////


  AllTrees(const std::string& input_path, const snapnum_type snapnum_first,
      const snapnum_type snapnum_last, const std::string& skipsnaps_filename)
      : subhalos_() {

    // Create subhalo objects
    get_subhalos(input_path, snapnum_first, snapnum_last, skipsnaps_filename);

  }


private:
  ///////////////////
  // PRIVATE TYPES //
  ///////////////////

  /** Internal type for subhalos. */
  struct internal_subhalo {

    // Some basic info
    sub_id_type id;
    snapnum_type snap;
    index_type index;
    index_type desc_index;
    index_type group_index;
    sub_len_type num_particles;
    real_type mass;
    real_type mass_history;
    uint8_t skip_snapshot;

    // Links to other subhalos.
    internal_subhalo* first_progenitor_ = nullptr;
    internal_subhalo* next_progenitor_ = nullptr;
    internal_subhalo* descendant_ = nullptr;
    internal_subhalo* first_subhalo_in_fof_group_ = nullptr;
    internal_subhalo* next_subhalo_in_fof_group_ = nullptr;
    internal_subhalo* last_progenitor_ = nullptr;
    internal_subhalo* main_leaf_progenitor_ = nullptr;
    internal_subhalo* root_descendant_ = nullptr;
    /** Constructor. */
    internal_subhalo(snapnum_type snap_, index_type index_,
        index_type desc_index_, index_type group_index_,
        sub_len_type num_particles_, real_type mass_, uint8_t skip_snapshot_)
        : id(-1), snap(snap_), index(index_), desc_index(desc_index_),
          group_index(group_index_), num_particles(num_particles_),
          mass(mass_), mass_history(0), skip_snapshot(skip_snapshot_),
          first_progenitor_(nullptr),
          next_progenitor_(nullptr),
          descendant_(nullptr),
          first_subhalo_in_fof_group_(nullptr),
          next_subhalo_in_fof_group_(nullptr),
          last_progenitor_(nullptr),
          main_leaf_progenitor_(nullptr),
          root_descendant_(nullptr) {
    }
    /** Default destructor. */
    ~internal_subhalo() = default;
  };

  //////////////////////////////
  // PRIVATE MEMBER FUNCTIONS //
  //////////////////////////////

  /** @brief Create subhalo objects and set all their member variables,
   *         except for their unique subhalo ID.
   */
  void get_subhalos(const std::string& input_path,
      const snapnum_type snapnum_first, const snapnum_type snapnum_last,
      const std::string& skipsnaps_filename) {

    // Check that subhalos_ structure has the correct size.
    assert(subhalos_.size() == 0);
    subhalos_.resize(snapnum_last+1);

    // Create list of valid snapshots.
    auto valid_snapnums = get_valid_snapnums(skipsnaps_filename,
        snapnum_first, snapnum_last);

    // ----------------- SNAPSHOT ITERATION 1

    // Initialize subhalo objects.
    std::cout << "Reading data and initializing subhalo objects...\n";
    WallClock wall_clock;
    for (auto snap_it = valid_snapnums.begin();
        snap_it != valid_snapnums.end(); ++snap_it) {
      auto cur_snapnum = *snap_it;

      // Create filename.
      std::stringstream tmp_stream;
      tmp_stream << input_path << "_" <<
          std::setfill('0') << std::setw(3) << cur_snapnum << ".hdf5";
      std::string desc_filename = tmp_stream.str();

      // Only proceed if current descendant file is non-empty.
      H5::H5File tmp_file(desc_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
      if (!H5Lexists(tmp_file.getId(), "/DescendantIndex", H5P_DEFAULT)) {
        tmp_file.close();
        continue;
      }
      tmp_file.close();

      // Read info from descendants file
      auto sub_len = read_dataset<sub_len_type>(desc_filename, "SubhaloLen");
      auto sub_mass = read_dataset<real_type>(desc_filename, "SubhaloMass");
      auto sub_grnr = read_dataset<index_type>(desc_filename, "SubhaloGrNr");
      auto desc_index = read_dataset<index_type>(desc_filename, "DescendantIndex");
      auto skip_snap = read_dataset<uint8_t>(desc_filename, "SkipSnapshot");

      // Make sure that current "row" of subhalos_ has the correct size
      uint32_t nsubs = sub_len.size();
      assert(subhalos_[cur_snapnum].size() == 0);
      subhalos_[cur_snapnum].resize(nsubs, nullptr);

      // Create new subhalo objects (containing at least 1 particle).
      for (uint32_t i = 0; i < nsubs; ++i) {
        if (sub_len[i] > 0) {
          subhalos_[cur_snapnum][i] = new internal_subhalo(cur_snapnum, i,
              desc_index[i], sub_grnr[i], sub_len[i], sub_mass[i], skip_snap[i]);
        }
      }
    }
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    // ----------------- SNAPSHOT ITERATION 2

    // Create progenitor/descendant links between subhalos.
    std::cout << "Creating progenitor/descendant links...\n";
    wall_clock.start();
    bool first_nonempty = true;
    for (auto snap_it = valid_snapnums.begin();
        snap_it+1 != valid_snapnums.end(); ++snap_it) {

      // For clarity:
      snapnum_type cur_snapnum = *snap_it;
      snapnum_type desc_snapnum_1 = *(snap_it+1);
      snapnum_type desc_snapnum_2 = -1;
      if (snap_it+2 != valid_snapnums.end()) {
        assert(snap_it+2 < valid_snapnums.end());  // test pointer arithmetic
        desc_snapnum_2 = *(snap_it+2);
      }

      // For the first non-empty snapshot, MassHistory equals Mass.
      auto& cur_snap = subhalos_[cur_snapnum];
      if (first_nonempty) {
        for (auto sub_it = cur_snap.begin(); sub_it != cur_snap.end(); ++sub_it) {
          auto cur_sub = *sub_it;
          if (cur_sub != nullptr)
            cur_sub->mass_history = cur_sub->mass;
        }
        first_nonempty = false;
      }

      // Create links between progenitors and descendants
      // and calculate MassHistory.
      for (auto sub_it = cur_snap.begin(); sub_it != cur_snap.end(); ++sub_it) {
        auto cur_prog = *sub_it;

        // Only proceed if current subhalo is valid and has a descendant
        if ((cur_prog == nullptr) || ((cur_prog->desc_index == -1)))
          continue;

        // Define link to descendant
        internal_subhalo* cur_desc;
        if (cur_prog->skip_snapshot == 0) {
          cur_desc = subhalos_[desc_snapnum_1][cur_prog->desc_index];
        }
        else {
          assert(cur_prog->skip_snapshot == 1);
          assert(desc_snapnum_2 != -1);
          cur_desc = subhalos_[desc_snapnum_2][cur_prog->desc_index];
        }
        cur_prog->descendant_ = cur_desc;

        // The progenitors of a subhalo are ordered by their mass history.
        // Put cur_sub in its proper "place."
        if (cur_prog->mass_history > cur_desc->mass_history) {
          // cur_prog is the largest so far; place at beginning of list
          if (cur_desc->first_progenitor_ != nullptr) {
            // The former first_progenitor becomes next_progenitor of cur_prog
            cur_prog->next_progenitor_ = cur_desc->first_progenitor_;
          }
          cur_desc->mass_history = cur_prog->mass_history;
          cur_desc->first_progenitor_ = cur_prog;
        }
        else {
          // Iterate over the next_progenitor link until we find
          // a progenitor with a smaller mass history.
          auto prev_prog = cur_desc->first_progenitor_;
          auto next_prog = cur_desc->first_progenitor_->next_progenitor_;
          while (true) {
            if (next_prog == nullptr) {
              // cur_prog is the smallest so far; add to end of list
              prev_prog->next_progenitor_ = cur_prog;
              break;
            }
            if (cur_prog->mass_history > next_prog->mass_history) {
              // Place cur_prog "between" prev_prog and cur_prog
              prev_prog->next_progenitor_ = cur_prog;
              cur_prog->next_progenitor_ = next_prog;
              break;
            }
            prev_prog = next_prog;
            next_prog = next_prog->next_progenitor_;
          }
        }
      }
      // Set final mass_history values for (first) descendant snapshot.
      auto& desc_snap = subhalos_[desc_snapnum_1];
      for (auto sub_it = desc_snap.begin(); sub_it != desc_snap.end(); ++sub_it) {
        auto cur_desc = *sub_it;
        if (cur_desc != nullptr)
          cur_desc->mass_history += cur_desc->mass;
      }
    }
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    // ----------------- SNAPSHOT ITERATION 3

    // Establish links between subhalos in the same FOF group, which are
    // also ordered by their mass history.
    std::cout << "Creating links within FoF groups...\n";
    wall_clock.start();
    for (auto snap_it = valid_snapnums.begin();
        snap_it != valid_snapnums.end(); ++snap_it) {

      // Keep main subhalo from each FoF group in this structure:
      std::vector<internal_subhalo*> main_sub_in_fof_group;

      // Iterate over subhalos
      auto cur_snapnum = *snap_it;
      auto& cur_snap = subhalos_[cur_snapnum];
      for (auto sub_it = cur_snap.begin(); sub_it != cur_snap.end(); ++sub_it) {
        auto cur_sub = *sub_it;
        if (cur_sub == nullptr)
          continue;

        // Resize structure if necessary
        if (static_cast<std::size_t>(cur_sub->group_index+1) >
            main_sub_in_fof_group.size())
          main_sub_in_fof_group.resize(cur_sub->group_index+1);

        auto& cur_main_sub = main_sub_in_fof_group[cur_sub->group_index];
        if (cur_main_sub == nullptr) {
          cur_main_sub = cur_sub;
          continue;
        }

        if (cur_sub->mass_history > cur_main_sub->mass_history) {
          // Link to former main subhalo in FoF group
          cur_sub->next_subhalo_in_fof_group_ = cur_main_sub;
          // Establish new main subhalo in FoF group
          cur_main_sub = cur_sub;
        }
        else {
          // Iterate over the next_subhalo_in_fof_group link until we find
          // a subhalo with a smaller mass history.
          auto prev_sub = cur_main_sub;
          auto next_sub = cur_main_sub->next_subhalo_in_fof_group_;
          while (true) {
            if (next_sub == nullptr) {
              // cur_sub is the smallest so far; add to end of list
              prev_sub->next_subhalo_in_fof_group_ = cur_sub;
              break;
            }
            if (cur_sub->mass_history > next_sub->mass_history) {
              // Place cur_prog "between" prev_prog and cur_prog
              prev_sub->next_subhalo_in_fof_group_ = cur_sub;
              cur_sub->next_subhalo_in_fof_group_ = next_sub;
              break;
            }
            prev_sub = next_sub;
            next_sub = next_sub->next_subhalo_in_fof_group_;
          }
        }
      }
      // Set first_subhalo_in_fof_group link
      for (auto sub_it = cur_snap.begin(); sub_it != cur_snap.end(); ++sub_it) {
        auto cur_sub = *sub_it;
        if (cur_sub != nullptr) {
          cur_sub->first_subhalo_in_fof_group_ =
              main_sub_in_fof_group[cur_sub->group_index];
        }
      }
    }
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  }


  //////////////////////////////
  // PRIVATE MEMBER VARIABLES //
  //////////////////////////////

  /** Auxiliary structure. For a given @a snapnum cur_snapnumsubfind_id,
   * returns a pointer to the corresponding internal_subhalo.
   */
  std::vector<std::vector<internal_subhalo*>> subhalos_;


};

