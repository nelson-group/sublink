#pragma once
/** @file ReadTreeHDF5.hpp
 * @brief Define a class for reading and traversing SubLink merger trees.
 *
 * See tree_test.cpp for an usage example.
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <cassert>

#include "GeneralHDF5.hpp"
#include "../Util/GeneralUtil.hpp"  // totally_ordered

/** @class Tree
 * @brief A class that loads a merger tree, representing it
 *        as a linked-list structure.
 */
class Tree {
private:
 // Predeclare internal type for subhalos.
 struct internal_subhalo;

public:
  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  /** @brief Type of this Tree. */
  typedef Tree tree_type;

  // Predeclaration of Subhalo class.
  class Subhalo;
  /** @brief Type of Subhalos. */
  typedef Subhalo subhalo_type;

  // Predeclaration of SnapshotIterator class.
  class SnapshotIterator;
  /** @brief Synonym for SnapshotIterator. */
  typedef SnapshotIterator snapshot_iterator;

  // Predeclaration of BranchIterator class.
  class BranchIterator;
  /** @brief Synonym for BranchIterator. */
  typedef BranchIterator branch_iterator;

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

  /** @brief Format of the subhalo data.
   *
   * Extra quantities can be chosen from the list found in
   * http://www.illustris-project.org/w/index.php/Merger_Trees#Tree_format.
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

  /** @brief Synonym for DataFormat. */
  typedef DataFormat data_format;

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  /** @brief Construct a Tree from an HDF5 file.
   * @param[in] treedir Directory containing the merger tree files.
   * @param[in] name Suffix of the merger tree filenames, e.g., "extended".
   * @param[in] filenum File number; -1 for the full, "concatenated"
   *                    merger tree file.
   */
  Tree(const std::string& treedir, const std::string& name, const int filenum)
      : subhalos_() {

    // Create filename for the given arguments.
    std::stringstream tmp_stream;
    if (filenum == -1)
      tmp_stream << treedir << "/" << name << ".hdf5";
    else
      tmp_stream << treedir << "/" << name << "." << filenum << ".hdf5";
    std::string treefilename = tmp_stream.str();

    // Create temporary subhalo map
    std::map<sub_id_type, internal_subhalo*> sub_map;
    create_subhalo_map(treefilename, sub_map);

    // Link subhalos
    std::cout << "Linking subhalos..." << std::endl;
    WallClock wall_clock;
    for (auto it = sub_map.begin(); it != sub_map.end(); ++it) {
      internal_subhalo* sub = it->second;

      // Increase size of subhalos_ as necessary, initializing
      // new elements to nullptr. Note that subhalos_[snapnum].size()
      // is not necessarily equal to the number of subhalos in the
      // corresponding Subfind catalog.
      if (static_cast<std::size_t>(sub->data_.SnapNum+1) > subhalos_.size())
        subhalos_.resize(sub->data_.SnapNum+1);
      if (static_cast<std::size_t>(sub->data_.SubfindID+1) >
          subhalos_[sub->data_.SnapNum].size())
        subhalos_[sub->data_.SnapNum].resize(sub->data_.SubfindID+1, nullptr);

      // Add subhalo to corresponding snapshot
      subhalos_[sub->data_.SnapNum][sub->data_.SubfindID] = sub;

      // Establish links between subhalos
      if (sub->data_.FirstProgenitorID != -1)
        sub->first_progenitor_ = sub_map[sub->data_.FirstProgenitorID];
      if (sub->data_.NextProgenitorID != -1)
        sub->next_progenitor_ = sub_map[sub->data_.NextProgenitorID];
      if (sub->data_.DescendantID != -1)
        sub->descendant_ = sub_map[sub->data_.DescendantID];
      if (sub->data_.FirstSubhaloInFOFGroupID != -1)
        sub->first_subhalo_in_fof_group_ = sub_map[sub->data_.FirstSubhaloInFOFGroupID];
      if (sub->data_.NextSubhaloInFOFGroupID != -1)
        sub->next_subhalo_in_fof_group_ = sub_map[sub->data_.NextSubhaloInFOFGroupID];
      if (sub->data_.LastProgenitorID != -1)
        sub->last_progenitor_ = sub_map[sub->data_.LastProgenitorID];
      if (sub->data_.MainLeafProgenitorID != -1)
        sub->main_leaf_progenitor_ = sub_map[sub->data_.MainLeafProgenitorID];
      if (sub->data_.RootDescendantID != -1)
        sub->root_descendant_ = sub_map[sub->data_.RootDescendantID];
    }
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  }

  /** Destructor. */
  ~Tree() {
    // Iterate over vector of vectors
    for (auto it1 = subhalos_.begin(); it1 != subhalos_.end(); ++it1) {
      for (auto it2 = (*it1).begin(); it2 != (*it1).end(); ++it2) {
        internal_subhalo* cur_sub = *it2;
        delete cur_sub;
      }
    }
  }

  /////////////
  // GENERAL //
  /////////////

  /** Return the number of snapshots in the tree.
   *
   * @note Equals @a snapnum_max+1, where @a snapnum_max is the maximum
   * snapshot number.
   */
  uint16_t num_snapshots() const {
    return subhalos_.size();
  }

  ///////////////
  // SNAPSHOTS //
  ///////////////

  /** @class Snapshot
   * @brief Class to represent a snapshot in the merger tree.
   */
  class Snapshot : private totally_ordered<Snapshot> {
  public:
    /** Construct invalid Snapshot. */
    Snapshot() : t_(nullptr), snap_(-1) {
    }

    /** Return snapshot number. */
    snapnum_type snapnum() const {
      return snap_;
    }

    /** Test whether this Snapshot and @a x are equal. */
    bool operator==(const Snapshot& x) const {
      return (t_ == x.t_) && (snap_ == x.snap_);
    }
    /** Test whether this Snapshot is less than @a x in the global order. */
    bool operator<(const Snapshot& x) const {
      // Compare t_, then snap_
      return std::tie(t_, snap_) < std::tie(x.t_, x.snap_);
    }

    /** Return an iterator pointing to the first subhalo in this snapshot. */
    SnapshotIterator begin() const {
      assert(is_valid());
      return SnapshotIterator(t_, snap_, 0);
    }
    /** Return iterator pointing to one-past-the-last subhalo in this
     * snapshot.
     */
    SnapshotIterator end() const {
      assert(is_valid());
      return SnapshotIterator(t_, snap_, t_->subhalos_[snap_].size());
    }

    /** Helper function to determine whether this Snapshot is valid. */
    bool is_valid() const {
      if (t_ != nullptr)
        return (snap_ >= 0) && (snap_ < t_->num_snapshots());
      return false;
    }

  private:
    friend class Tree;
    // Pointer back to the Tree containing this Snapshot.
    Tree* t_;
    // Snapshot number.
    snapnum_type snap_;
    /** Private constructor. */
    Snapshot(const tree_type* t, snapnum_type snap)
        : t_(const_cast<tree_type*>(t)), snap_(snap) {
    }
  };

  /** @brief Return the Snapshot object given by @a snapnum.
   * @pre 0 <= @a snapnum < num_snapshots()
   */
  Snapshot snapshot(snapnum_type snapnum) const {
    assert((snapnum >= 0) && (snapnum < num_snapshots()));
    return Snapshot(this, snapnum);
  }

  //////////////
  // SUBHALOS //
  //////////////

  /** @class Subhalo
   * @brief Class to represent subhalos in the merger tree.
   */
  class Subhalo : private totally_ordered<Subhalo> {
  public:
    /** Construct an invalid Subhalo. */
    Subhalo() : t_(nullptr), snap_(-1), idx_(-1) {
    }

    /** Return the snapshot number of this Subhalo. */
    snapnum_type snapnum() const {
      return snap_;
    }
    /** Return the Subfind ID of this Subhalo. */
    index_type index() const {
      return idx_;
    }

    /** Test whether this Subhalo and @a x are equal. */
    bool operator==(const Subhalo& x) const {
      return (t_ == x.t_) && (snap_ == x.snap_) && (idx_ == x.idx_);
    }
    /** Test whether this Subhalo is less than @a x in the global order.
     * @note An alternative ordering is given by SubhaloID, which orders
     *       subhalos in a depth-first fashion in the merger tree.
     */
    bool operator<(const Subhalo& x) const {
      // Compare t_, then snap_, then idx_
      return std::tie(t_, snap_, idx_) < std::tie(x.t_, x.snap_, x.idx_);
    }

    /** Return reference to data as found in the merger tree. */
    const DataFormat& data() const {
      return fetch()->data_;
    }

    // The member functions below return "linked" Subhalo objects.
    Subhalo first_progenitor() const {
      return ptr_to_sub(fetch()->first_progenitor_);
    }
    Subhalo next_progenitor() const {
      return ptr_to_sub(fetch()->next_progenitor_);
    }
    Subhalo descendant() const {
      return ptr_to_sub(fetch()->descendant_);
    }
    Subhalo first_subhalo_in_fof_group() const {
      return ptr_to_sub(fetch()->first_subhalo_in_fof_group_);
    }
    Subhalo next_subhalo_in_fof_group() const {
      return ptr_to_sub(fetch()->next_subhalo_in_fof_group_);
    }
    Subhalo last_progenitor() const {
      return ptr_to_sub(fetch()->last_progenitor_);
    }
    Subhalo main_leaf_progenitor() const {
      return ptr_to_sub(fetch()->main_leaf_progenitor_);
    }
    Subhalo root_descendant() const {
      return ptr_to_sub(fetch()->root_descendant_);
    }

    /** Return a BranchIterator pointing to the current Subhalo. */
    BranchIterator branch_begin() const {
      auto p = fetch();
      return BranchIterator(t_, p, p->main_leaf_progenitor_);
    }
    /** Return a BranchIterator pointing to "one past" the main leaf
     * progenitor of this Subhalo.
     */
    BranchIterator branch_end() const {
      auto p = fetch();
      return BranchIterator(t_, nullptr, p->main_leaf_progenitor_);
    }

    /** Helper function to determine whether this Subhalo is valid. */
    bool is_valid() const {
      return (t_ != nullptr) &&
          (snap_ >= 0) && (static_cast<std::size_t>(snap_) < t_->num_snapshots()) &&
          (idx_ >= 0)  && (static_cast<std::size_t>(idx_)  < t_->subhalos_[snap_].size()) &&
          (t_->subhalos_[snap_][idx_] != nullptr);
    }

  private:
    friend class Tree;
    // Pointer back to the Tree containing this Subhalo.
    Tree* t_;
    // Snapshot number where this subhalo is found.
    snapnum_type snap_;
    // Index of this subhalo in the Subfind catalog.
    index_type idx_;
    /** Private constructor. */
    Subhalo(const tree_type* t, snapnum_type snap, index_type idx)
        : t_(const_cast<tree_type*>(t)), snap_(snap), idx_(idx) {
    }
    /** Return a pointer to the corresponding internal_subhalo. */
    internal_subhalo* fetch() const {
      assert(is_valid());
      return t_->subhalos_[snap_][idx_];
    }
    /** Return a Subhalo object for a given internal_subhalo*.
     * @param[in] p Pointer to @a internal_subhalo object.
     * @pre @a t_ != @a nullptr.
     */
    Subhalo ptr_to_sub(internal_subhalo* p) const {
      if (p != nullptr)
        return Subhalo(t_, p->data_.SnapNum, p->data_.SubfindID);
      return Subhalo();
    }
  };

  /** @brief Return the Subhalo object given by @a snapnum and @a subfind_id.
   * @pre 0 <= @a snapnum < num_snapshots() && @a idx >= 0
   *
   * If the requested subhalo does not exist in the tree, an invalid Subhalo
   * object is returned.
   */
  Subhalo subhalo(snapnum_type snapnum, index_type idx) const {
    assert((snapnum >= 0) && (snapnum < num_snapshots()) && (idx >= 0));
    if ((static_cast<std::size_t>(idx) >= subhalos_[snapnum].size()) ||
        (subhalos_[snapnum][idx] == nullptr))
      return Subhalo();
    return Subhalo(this, snapnum, idx);
  }

  ///////////////
  // ITERATORS //
  ///////////////

  /** @class Tree::SnapshotIterator
   * @brief Iterator class for Subhalos within the same snapshot.
   *
   * @note Only iterates over subhalos that exist in the merger tree.
   */
  class SnapshotIterator : private totally_ordered<SnapshotIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Subhalo value_type;
    /** Type of pointers to elements. */
    typedef Subhalo* pointer;
    /** Type of references to elements. */
    typedef Subhalo& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid SnapshotIterator. */
    SnapshotIterator() : t_(nullptr), snap_(-1), idx_(-1) {
    }

    /** Method to dereference a SnapshotIterator.
     * @pre Subhalo must be valid.
     */
    Subhalo operator*() const {
      assert(is_valid());
      return Subhalo(t_, snap_, idx_);
    }
    /** Method to increment a SnapshotIterator. */
    SnapshotIterator& operator++() {
      ++idx_;
      fix();
      return *this;
    }
    /** Method to compare two SnapshotIterators.
     * @note Invalid SnapshotIterators are always equal to each other.
     */
    bool operator==(const SnapshotIterator& x) const {
      return (t_ == x.t_) && (snap_ == x.snap_) && (idx_ == x.idx_);
    }

   private:
    friend class Tree;
    Tree* t_;
    snapnum_type snap_;
    index_type idx_;
    /** Private constructor. */
    SnapshotIterator(const tree_type* t, snapnum_type snap, index_type idx)
        : t_(const_cast<tree_type*>(t)), snap_(snap), idx_(idx) {
      fix();
    }
    /** Helper function to determine whether this iterator is valid. */
    bool is_valid() const {
      return (t_ != nullptr) &&
          (snap_ >= 0) && (static_cast<std::size_t>(snap_) < t_->num_snapshots()) &&
          (idx_ >= 0)  && (static_cast<std::size_t>(idx_)  < t_->subhalos_[snap_].size()) &&
          (t_->subhalos_[snap_][idx_] != nullptr);
    }
    /** Function to advance an iterator to a valid position or end(). */
    void fix() {
      assert(t_ != nullptr);
      assert((snap_ >= 0) && (static_cast<std::size_t>(snap_) < t_->num_snapshots()));
      assert(idx_ >= 0);
      auto nsubs = t_->subhalos_[snap_].size();
      // If idx_ >= nsubs, return end().
      if (static_cast<std::size_t>(idx_) >= nsubs) {
        idx_ = nsubs;
        return;
      }
      // Advance until a valid Subhalo or end() is reached.
      while ((t_->subhalos_[snap_][idx_] == nullptr) &&
          (static_cast<std::size_t>(idx_) < nsubs))
        ++idx_;
    }
  };

  /** @class Tree::BranchIterator
   * @brief Iterator for Subhalos along the main branch of a given Subhalo.
   *
   * RI(*this): cur_sub_->main_leaf_progenitor_ == main_leaf_
   *
   * Implementation note: the main_leaf_ variable is redundant,
   * but is included for consistency checks.
   */
  class BranchIterator : private totally_ordered<BranchIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Subhalo value_type;
    /** Type of pointers to elements. */
    typedef Subhalo* pointer;
    /** Type of references to elements. */
    typedef Subhalo& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid BranchIterator. */
    BranchIterator() : t_(nullptr), cur_sub_(nullptr), main_leaf_(nullptr) {
    }

    /** Method to dereference a BranchIterator.
     * @pre Iterator must be valid.
     */
    Subhalo operator*() const {
      assert(is_valid());
      return Subhalo(t_, cur_sub_->data_.SnapNum, cur_sub_->data_.SubfindID);
    }
    /** Method to increment a BranchIterator.
     * @pre Iterator must be valid.
     */
    BranchIterator& operator++() {
      assert(is_valid());
      cur_sub_ = cur_sub_->first_progenitor_;  // can be nullptr
      return *this;
    }
    /** Method to compare two BranchIterators.
     * @note Invalid BranchIterators are always equal to each other.
     */
    bool operator==(const BranchIterator& x) const {
      return (t_ == x.t_) && (main_leaf_ == x.main_leaf_) &&
          (cur_sub_ == x.cur_sub_);
    }

   private:
    friend class Tree;
    // Pointer back to the Tree container.
    Tree* t_;
    // Pointer to the current subhalo.
    internal_subhalo* cur_sub_;
    // Pointer to the main leaf progenitor (only used for consistency checks).
    internal_subhalo* main_leaf_;
    /** Private constructor. */
    BranchIterator(const tree_type* t, internal_subhalo* cur_sub,
        internal_subhalo* main_leaf)
        : t_(const_cast<tree_type*>(t)), cur_sub_(cur_sub),
          main_leaf_(main_leaf) {
    }
    /** Helper function to determine whether this iterator is valid. */
    bool is_valid() const {
      return (t_ != nullptr) && (main_leaf_ != nullptr) && (cur_sub_ != nullptr) &&
          (cur_sub_->main_leaf_progenitor_ == main_leaf_);
    }
  };

private:
  ///////////////////
  // PRIVATE TYPES //
  ///////////////////

  /** Internal type for subhalos. */
  struct internal_subhalo {
    // Fields from "minimal" data format.
    DataFormat data_;
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
    internal_subhalo(const DataFormat data)
        : data_(data),
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

  /** @brief Read data from HDF5 file and return a temporary mapping between
   *         SubhaloIDs and pointers to internal_subhalos.
   */
  template <typename MAP>
  void create_subhalo_map(const std::string& treefilename, MAP& sub_map) {
    // Read data
    std::cout << "Reading data from file " << treefilename << "\n";
    WallClock wall_clock;
    auto SubhaloID = read_dataset<sub_id_type>(treefilename, "SubhaloID");
    auto SubhaloIDRaw = read_dataset<sub_id_type>(treefilename, "SubhaloIDRaw");
    auto LastProgenitorID = read_dataset<sub_id_type>(treefilename, "LastProgenitorID");
    auto MainLeafProgenitorID = read_dataset<sub_id_type>(treefilename, "MainLeafProgenitorID");
    auto RootDescendantID = read_dataset<sub_id_type>(treefilename, "RootDescendantID");
    auto TreeID = read_dataset<tree_id_type>(treefilename, "TreeID");
    auto SnapNum = read_dataset<snapnum_type>(treefilename, "SnapNum");
    auto FirstProgenitorID = read_dataset<sub_id_type>(treefilename, "FirstProgenitorID");
    auto NextProgenitorID = read_dataset<sub_id_type>(treefilename, "NextProgenitorID");
    auto DescendantID = read_dataset<sub_id_type>(treefilename, "DescendantID");
    auto FirstSubhaloInFOFGroupID = read_dataset<sub_id_type>(treefilename, "FirstSubhaloInFOFGroupID");
    auto NextSubhaloInFOFGroupID = read_dataset<sub_id_type>(treefilename, "NextSubhaloInFOFGroupID");
    auto NumParticles = read_dataset<uint32_t>(treefilename, "NumParticles");
    auto Mass = read_dataset<real_type>(treefilename, "Mass");
    auto MassHistory = read_dataset<real_type>(treefilename, "MassHistory");
    auto SubfindID = read_dataset<index_type>(treefilename, "SubfindID");
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    // Create internal_subhalo objects, storing pointers to them in a map.
    wall_clock.start();
    std::cout << "Creating subhalo map...\n";
    uint64_t nrows = SubhaloID.size();
    for (uint64_t rownum = 0; rownum < nrows; ++rownum) {
      sub_map.emplace(SubhaloID[rownum], new internal_subhalo(
          DataFormat(
              SubhaloID[rownum],
              SubhaloIDRaw[rownum],
              LastProgenitorID[rownum],
              MainLeafProgenitorID[rownum],
              RootDescendantID[rownum],
              TreeID[rownum],
              SnapNum[rownum],
              FirstProgenitorID[rownum],
              NextProgenitorID[rownum],
              DescendantID[rownum],
              FirstSubhaloInFOFGroupID[rownum],
              NextSubhaloInFOFGroupID[rownum],
              NumParticles[rownum],
              Mass[rownum],
              MassHistory[rownum],
              SubfindID[rownum])));
    }
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  }

  //////////////////////////////
  // PRIVATE MEMBER VARIABLES //
  //////////////////////////////

  /** Auxiliary structure. For a given @a snapnum and @a subfind_id,
   * returns a pointer to the corresponding internal_subhalo.
   */
  std::vector<std::vector<internal_subhalo*>> subhalos_;
};
