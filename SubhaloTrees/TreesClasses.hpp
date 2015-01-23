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
#include <algorithm>  // sort

#include "../InputOutput/GeneralHDF5.hpp"
#include "../Util/SnapshotUtil.hpp"
#include "../Util/GeneralUtil.hpp"
#include "../Util/TreeTypes.hpp"

static constexpr int64_t pow_10_8  = 100000000;
static constexpr int64_t pow_10_12 = 1000000000000;
static constexpr int64_t pow_10_16 = 10000000000000000;

/** @class AllTrees
 * @brief A class for constructing merger trees.
 */
class AllTrees {
private:
  // Predeclare internal types.
  struct internal_subhalo;
  struct internal_tree;

public:
  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  /** @brief Format of the subhalo data, consisting of the "minimal" fields. */
  struct DataFormat {
    // Fields corresponding to "minimal" data format.
    sub_id_type SubhaloID;
    sub_id_type SubhaloIDRaw;
    sub_id_type LastProgenitorID;
    sub_id_type MainLeafProgenitorID;
    sub_id_type RootDescendantID;
    tree_id_type TreeID;
    // Although 2 bytes should be enough, this entry occupies 8 bytes anyway:
    int64_t SnapNum;
    sub_id_type FirstProgenitorID;
    sub_id_type NextProgenitorID;
    sub_id_type DescendantID;
    sub_id_type FirstSubhaloInFOFGroupID;
    sub_id_type NextSubhaloInFOFGroupID;
    sub_len_type NumParticles;
    real_type Mass;
    real_type MassHistory;
    index_type SubfindID;

    // Forbid default construction.
    DataFormat() = delete;

    /** @brief Construct a DataFormat object from a given @a internal_subhalo.
     * @param[in] sub A pointer to the @a internal_subhalo object.
     * @pre @a sub != @a nullptr.
     */
    DataFormat(const internal_subhalo * const & sub)
        : SubhaloID(sub->id),
          SubhaloIDRaw(pow_10_12*sub->snap + sub->index),
          LastProgenitorID(sub->last_progenitor_ == nullptr ? -1 : sub->last_progenitor_->id),
          MainLeafProgenitorID(sub->main_leaf_progenitor_ == nullptr ? -1 : sub->main_leaf_progenitor_->id),
          RootDescendantID(sub->root_descendant_ == nullptr ? -1 : sub->root_descendant_->id),
          TreeID(sub->tree_ == nullptr ? -1 : sub->tree_->id),
          SnapNum(sub->snap),
          FirstProgenitorID(sub->first_progenitor_ == nullptr ? -1 : sub->first_progenitor_->id),
          NextProgenitorID(sub->next_progenitor_ == nullptr ? -1 : sub->next_progenitor_->id),
          DescendantID(sub->descendant_ == nullptr ? -1 : sub->descendant_->id),
          FirstSubhaloInFOFGroupID(sub->first_subhalo_in_fof_group_ == nullptr ? -1 : sub->first_subhalo_in_fof_group_->id),
          NextSubhaloInFOFGroupID(sub->next_subhalo_in_fof_group_ == nullptr ? -1 : sub->next_subhalo_in_fof_group_->id),
          NumParticles(sub->num_particles),
          Mass(sub->mass),
          MassHistory(sub->mass_history),
          SubfindID(sub->index) {
    }
  };

  /** @brief Format of the subhalo data, HDF5 version. */
  static H5::CompType H5DataFormat() {
    H5::CompType mtype(sizeof(DataFormat));
    mtype.insertMember("SubhaloID", HOFFSET(DataFormat, SubhaloID), H5::PredType::NATIVE_INT64);
    mtype.insertMember("SubhaloIDRaw", HOFFSET(DataFormat, SubhaloIDRaw), H5::PredType::NATIVE_INT64);
    mtype.insertMember("LastProgenitorID", HOFFSET(DataFormat, LastProgenitorID), H5::PredType::NATIVE_INT64);
    mtype.insertMember("MainLeafProgenitorID", HOFFSET(DataFormat, MainLeafProgenitorID), H5::PredType::NATIVE_INT64);
    mtype.insertMember("RootDescendantID", HOFFSET(DataFormat, RootDescendantID), H5::PredType::NATIVE_INT64);
    mtype.insertMember("TreeID", HOFFSET(DataFormat, TreeID), H5::PredType::NATIVE_INT64);
    // Although 2 bytes should be enough, this entry occupies 8 bytes anyway:
    mtype.insertMember("SnapNum", HOFFSET(DataFormat, SnapNum), H5::PredType::NATIVE_INT64);
    mtype.insertMember("FirstProgenitorID", HOFFSET(DataFormat, FirstProgenitorID), H5::PredType::NATIVE_INT64);
    mtype.insertMember("NextProgenitorID", HOFFSET(DataFormat, NextProgenitorID), H5::PredType::NATIVE_INT64);
    mtype.insertMember("DescendantID", HOFFSET(DataFormat, DescendantID), H5::PredType::NATIVE_INT64);
    mtype.insertMember("FirstSubhaloInFOFGroupID", HOFFSET(DataFormat, FirstSubhaloInFOFGroupID), H5::PredType::NATIVE_INT64);
    mtype.insertMember("NextSubhaloInFOFGroupID", HOFFSET(DataFormat, NextSubhaloInFOFGroupID), H5::PredType::NATIVE_INT64);
    mtype.insertMember("NumParticles", HOFFSET(DataFormat, NumParticles), H5::PredType::NATIVE_UINT32);
    mtype.insertMember("Mass", HOFFSET(DataFormat, Mass), H5::PredType::NATIVE_FLOAT);
    mtype.insertMember("MassHistory", HOFFSET(DataFormat, MassHistory), H5::PredType::NATIVE_FLOAT);
    mtype.insertMember("SubfindID", HOFFSET(DataFormat, SubfindID), H5::PredType::NATIVE_INT32);
    return mtype;
  }

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  /** Constructor. Does all the work. */
  AllTrees(const std::string& input_path, const std::string& output_path,
      const snapnum_type snapnum_first, const snapnum_type snapnum_last,
      const std::string& skipsnaps_filename)
      : subhalos_(), trees_() {

    // Create subhalo objects.
    get_subhalos(input_path, snapnum_first, snapnum_last, skipsnaps_filename);

    // Construct merger trees.
    first_pass();
    second_pass();

    // Assign IDs and write to files.
    write_to_files(output_path);
  }

  /** Destructor. */
  ~AllTrees() {
    for (auto it = trees_.begin(); it != trees_.end(); ++it)
      delete *it;

    for (auto it1 = subhalos_.begin(); it1 != subhalos_.end(); ++it1)
      for (auto it2 = (*it1).begin(); it2 != (*it1).end(); ++it2)
        delete *it2;
  }

private:
  ///////////////////
  // PRIVATE TYPES //
  ///////////////////

  /** @brief Internal type for subhalos. */
  struct internal_subhalo {
    // Some basic subhalo info
    sub_id_type id;
    snapnum_type snap;
    index_type index;
    index_type desc_index;
    index_type group_index;
    sub_len_type num_particles;
    real_type mass;
    // A double here, but ultimately converted to float:
    double mass_history;
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

    // Link to tree containing this subhalo.
    internal_tree* tree_ = nullptr;

    /** Forbid default constructor. */
    internal_subhalo() = delete;

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
          root_descendant_(nullptr),
          tree_(nullptr) {
    }

    /** Default destructor. */
    ~internal_subhalo() = default;
  };

  /** @brief Internal type for trees. */
  struct internal_tree {
    // Unique ID of this tree.
    tree_id_type id;
    // Subhalos belonging to this tree.
    std::vector<internal_subhalo*> subhalos;
    /** Default constructor. */
    internal_tree() : id(-1), subhalos() {
    }

    /** @brief Add a subhalo and all its progenitors to the tree
     *         in a depth-first fashion.
     *
     * @param[in] sub Pointer to an internal_subhalo object, such that
     *                all of its progenitors are added to *this.
     * @pre @a sub != nullptr
     * @note Also sets the @a root_descendant, @a last_progenitor, and
     * @a main_leaf_progenitor links along the way.
     *
     * @note Empty subhalos (i.e., with zero particles of the desired type)
     * cannot have "genealogic" links (i.e., progenitor or descendant).
     * Therefore, we cannot "bump" into such subhalos by recursively using
     * this function.
     */
    void add_recursive(internal_subhalo* sub) {
      assert(sub != nullptr);

      // Add current subhalo to current tree.
      sub->tree_ = this;
      subhalos.push_back(sub);

      // Assign root descendant.
      auto desc = sub->descendant_;
      if (desc == nullptr)
        sub->root_descendant_ = sub;
      else {
        // Root descendant is the same as for the descendant
        assert(desc->root_descendant_ != nullptr);
        sub->root_descendant_ = desc->root_descendant_;
      }

      // Add all progenitors to tree.
      auto first_prog = sub->first_progenitor_;
      if (first_prog == nullptr)
        sub->main_leaf_progenitor_ = sub;
      else {
        // Add first progenitor to tree, along with all its progenitors.
        add_recursive(first_prog);

        // Main leaf progenitor is the same as for the main progenitor.
        assert(first_prog->main_leaf_progenitor_ != nullptr);
        sub->main_leaf_progenitor_ = first_prog->main_leaf_progenitor_;

        // Add the other progenitors recursively.
        for (auto next_prog = first_prog->next_progenitor_;
             next_prog != nullptr; next_prog = next_prog->next_progenitor_)
          add_recursive(next_prog);
      }

      // The last progenitor of the current subhalo is the last object
      // added to the tree.
      sub->last_progenitor_ = subhalos.back();
    }
  };

  //////////////////////////////
  // PRIVATE MEMBER FUNCTIONS //
  //////////////////////////////

  /** @brief Create subhalo objects and set all their member variables,
   *         except for their unique subhalo ID and their root_descendant,
   *         last_progenitor and main_leaf_progenitor pointers.
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

      // Create new subhalo objects.
      for (uint32_t i = 0; i < nsubs; ++i) {
        // Only if subhalo contains at least 1 particle.
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
          main_sub_in_fof_group.resize(cur_sub->group_index+1, nullptr);

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

  /** @brief Merge two subtrees.
   *
   * Add subhalos from the smaller tree to the larger tree.
   * Then empty the smaller tree, which is ignored and deleted
   * while writing to files.
   */
  void merge_trees(internal_tree* tree1, internal_tree* tree2) {

    // Make sure that tree1 is the largest one.
    if (tree2->subhalos.size() > tree1->subhalos.size()) {
      auto tree_aux = tree1;
      tree1 = tree2;
      tree2 = tree_aux;
      assert(tree1->subhalos.size() > tree2->subhalos.size());
    }

    // Transfer subhalos from tree2 to tree1.
    for (auto sub_it = tree2->subhalos.begin();
         sub_it != tree2->subhalos.end(); ++sub_it) {
      auto cur_sub = *sub_it;
      cur_sub->tree_ = tree1;
      tree1->subhalos.push_back(cur_sub);
    }

    // Empty smaller tree
    std::vector<internal_subhalo*>().swap(tree2->subhalos);
    assert(tree2->subhalos.size() == 0);
  }

  /** @brief Construct merger trees. First stage.
   *
   * As a first approximation, each tree is formed by the subtrees
   * which are rooted in subhalos that belong to the same FoF group
   * at the latest snapshot.
   */
  void first_pass() {
    WallClock wall_clock;
    std::cout << "Constructing trees. First pass...\n";

    // Store tree corresponding to each FoF group in this vector.
    std::vector<internal_tree*> tmp_trees;

    // Iterate over subhalos from last snapshot.
    auto& last_snap = subhalos_.back();
    for (auto sub_it = last_snap.begin(); sub_it != last_snap.end(); ++sub_it) {
      auto cur_sub = *sub_it;

      // Only proceed if subhalo is valid.
      if (cur_sub == nullptr)
        continue;
      assert(cur_sub->num_particles > 0);

      // Resize vector if necessary
      assert(cur_sub->group_index >= 0);
      uint32_t group_index = cur_sub->group_index;
      if (group_index+1 > tmp_trees.size())
        tmp_trees.resize(group_index+1, nullptr);

      // Add contribution from subhalo to corresponding tree.
      if (tmp_trees[group_index] == nullptr)
        tmp_trees[group_index] = new internal_tree();
      tmp_trees[group_index]->add_recursive(cur_sub);
    }

    // Add trees to private data structure.
    for (auto tree_it = tmp_trees.begin(); tree_it != tmp_trees.end();
        ++tree_it) {
      auto cur_tree = *tree_it;
      if (cur_tree != nullptr)
        trees_.push_back(cur_tree);
    }
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  }

  /** @brief Construct merger trees. Second stage.
   *
   * Not all subhalos are "connected" to a subhalo from the last snapshot,
   * so we deal with those cases here. For each subhalo, we make the
   * following comparison with the main subhalo from its FoF group:
   *
   * 1) If one or both subhalos do not belong to a tree yet, we build
   * the corresponding subtrees and make sure that they all end up in the
   * same tree.
   *
   * 2) If the two subhalos belong to different trees, we merge them.
   *
   */
  void second_pass() {
    WallClock wall_clock;
    std::cout << "Constructing trees. Second pass...\n";

    // Iterate over snapshots in reverse.
    for (auto rit = subhalos_.rbegin(); rit != subhalos_.rend(); ++rit) {
      auto& cur_snap = *rit;
      // Iterate over subhalos in snapshot.
      for (auto sub_it = cur_snap.begin(); sub_it != cur_snap.end(); ++sub_it) {
        auto cur_sub = *sub_it;
        // Only proceed if subhalo is valid.
        if (cur_sub == nullptr)
          continue;
        assert(cur_sub->num_particles > 0);
        assert(cur_sub->first_subhalo_in_fof_group_ != nullptr);
        // Check if the main subhalo in the FoF group (which can be
        // cur_sub itself) already belongs to a tree.
        // If not, create a new tree rooted in the main subhalo.
        if (cur_sub->first_subhalo_in_fof_group_->tree_ == nullptr) {
          auto new_tree = new internal_tree();
          new_tree->add_recursive(cur_sub->first_subhalo_in_fof_group_);
          trees_.push_back(new_tree);
        }
        // Check if current subhalo already belongs to a tree.
        // If not, just add to tree of main and continue.
        if (cur_sub->tree_ == nullptr) {
          cur_sub->first_subhalo_in_fof_group_->tree_->add_recursive(cur_sub);
          continue;
        }
        // If cur_sub and the main subhalo belong to different trees,
        // we merge them.
        if (cur_sub->first_subhalo_in_fof_group_->tree_ != cur_sub->tree_)
          merge_trees(cur_sub->first_subhalo_in_fof_group_->tree_,
                      cur_sub->tree_);
      }
    }
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  }

  /** @brief Assign unique IDs and write to files. */
  void write_to_files(const std::string& writepath) {
    // Sort trees by decreasing size.
    std::cout << "Sorting trees by size...\n";
    WallClock wall_clock;
    std::sort(trees_.begin(), trees_.end(),
        [&](const internal_tree * const & tree1,
            const internal_tree * const & tree2) {
                return tree1->subhalos.size() > tree2->subhalos.size();
    });
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    // Remove (and delete) empty trees.
    while (trees_.back()->subhalos.size() == 0) {
      delete trees_.back();
      trees_.pop_back();
    }
    uint64_t ntrees = trees_.size();

    // Size of first tree determines approximate file size.
    // We make sure that no tree file is larger than the first one.
    // (alternatively, set max_nsubs_per_file to desired value)
    uint64_t max_nsubs_per_file = trees_.front()->subhalos.size();

    // Assign tree and subhalo IDs and determine which are the first
    // and last trees that go into each file.
    std::cout << "Assigning unique IDs...\n";
    wall_clock.start();
    std::vector<uint64_t> last_tree_in_file;
    std::vector<uint64_t> nsubs_per_file;
    int16_t filenum = 0;
    uint64_t tree_count = 0;
    uint64_t nsubs_in_cur_file = 0;
    for (uint64_t tree_index = 0; tree_index < ntrees; ++tree_index) {
      auto cur_tree = trees_[tree_index];
      // Assign IDs
      cur_tree->id = filenum*pow_10_16 + tree_count*pow_10_8;
      uint64_t subhalo_count = 0;
      for (auto sub_it = cur_tree->subhalos.begin();
          sub_it != cur_tree->subhalos.end(); ++sub_it) {
        auto cur_sub = *sub_it;
        assert(cur_sub != nullptr);
        cur_sub->id = cur_tree->id + subhalo_count;
        ++subhalo_count;
      }
      ++tree_count;
      // Check if this is the last tree that goes into the current file.
      nsubs_in_cur_file += cur_tree->subhalos.size();
      if (nsubs_in_cur_file + cur_tree->subhalos.size() >= max_nsubs_per_file) {
        last_tree_in_file.push_back(tree_index);
        nsubs_per_file.push_back(nsubs_in_cur_file);
        nsubs_in_cur_file = 0;
        tree_count = 0;
        ++filenum;
      }
    }
    // Put remaining trees, if any, into one last file.
    if (last_tree_in_file.back() != ntrees-1) {
      last_tree_in_file.push_back(ntrees-1);
      nsubs_per_file.push_back(nsubs_in_cur_file);
    }
    // Determine the first tree that goes into each file.
    std::vector<uint64_t> first_tree_in_file = {0};
    uint16_t nfiles = last_tree_in_file.size();
    for (filenum = 0; filenum < nfiles-1; ++filenum) {
      first_tree_in_file.push_back(last_tree_in_file[filenum]+1);
    }
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    // Write to files.
    std::cout << "Writing " << ntrees << " trees to " << nfiles << " files...\n";
    wall_clock.start();
    for (filenum = 0; filenum < nfiles; ++filenum) {

      // Create filename
      std::stringstream tmp_stream;
      tmp_stream << writepath << "." << filenum << ".hdf5";
      std::string writefilename = tmp_stream.str();

      // Create array with tree data.
      std::vector<DataFormat> treedata;
      treedata.reserve(nsubs_per_file[filenum]);
      for (uint64_t tree_index = first_tree_in_file[filenum];
           tree_index <= last_tree_in_file[filenum]; ++tree_index) {
        assert(tree_index < trees_.size());
        auto& cur_tree = trees_[tree_index];
        for (auto sub_it = cur_tree->subhalos.begin();
            sub_it != cur_tree->subhalos.end(); ++sub_it) {
          auto cur_sub = *sub_it;
          assert(cur_sub != nullptr);
          treedata.emplace_back(cur_sub);
        }
      }
      assert(treedata.size() == nsubs_per_file[filenum]);

      // Write to file.
      H5::H5File writefile(writefilename, H5F_ACC_TRUNC);
      add_array(writefile, treedata, "Tree", H5DataFormat());
      writefile.close();
    }

    std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  }

  //////////////////////////////
  // PRIVATE MEMBER VARIABLES //
  //////////////////////////////

  /** For a given @a snapnum and @a subfind_id,
   * @a subhalos_[snapnum][subfind_id] returns
   * a pointer to the corresponding internal_subhalo. */
  std::vector<std::vector<internal_subhalo*>> subhalos_;

  /** A vector with all the trees. */
  std::vector<internal_tree*> trees_;
};
