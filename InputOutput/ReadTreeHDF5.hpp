#pragma once
/** @file ReadTreeHDF5.hpp
 * @brief Define a class for reading merger trees.
 */

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>

#include "ReadWriteHDF5.hpp"

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

/** @class Tree
 * @brief A class used to traverse and extract information from
 *        SubLink merger trees.
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

  // Define data format.
  struct data_format {
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
    data_format(sub_id_type SubhaloID_,
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

  /** Construct a Tree from an HDF5 file. */
  Tree(const std::string& treedir, const std::string& name, const int filenum)
      : subhalos_() {

    // Create filename
    std::stringstream tmp_stream;
    if (filenum == -1)
      tmp_stream << treedir << "/" << name << ".hdf5";
    else
      tmp_stream << treedir << "/" << name << "." << filenum << ".hdf5";
    std::string treefilename = tmp_stream.str();

    // Create subhalo map
    std::map<sub_id_type, internal_subhalo*> sub_map;
    create_subhalo_map(treefilename, sub_map);

    // Link subhalos
    std::cout << "Linking subhalos..." << std::endl;
    for (auto it = sub_map.begin(); it != sub_map.end(); ++it) {
      internal_subhalo* sub = it->second;

      // Increase size of subhalos_ if necessary.
      if (static_cast<std::size_t>(sub->data_.SnapNum+1) > subhalos_.size())
        subhalos_.resize(sub->data_.SnapNum+1);
      if (static_cast<std::size_t>(sub->data_.SubfindID+1) >
          subhalos_[sub->data_.SnapNum].size())
        subhalos_[sub->data_.SnapNum].resize(sub->data_.SubfindID+1);

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
  // General //
  /////////////

  // num_snapshots, num_subhalos, etc.

  //////////////
  // SUBHALOS //
  //////////////

  /** @class Subhalo
   * @brief Class to represent subhalos in the merger tree.
   */
  class Subhalo {
  public:

    /** Construct an invalid Subhalo. */
    Subhalo() : t_(nullptr), snap_(-1), idx_(-1) {
    }

  private:
    friend class Tree;
    // Pointer back to the Tree containing this Subhalo.
    Tree* t_;
    // Snapshot number where this subhalo is found.
    snapnum_type snap_;
    // Index of this subhalo in the Subfind catalog.
    index_type idx_;
  };

private:
  ///////////////////
  // PRIVATE TYPES //
  ///////////////////

  /** Internal type for subhalos. */
  struct internal_subhalo {
    // Fields from "minimal" data format.
    data_format data_;

    // Links to other subhalos.
    internal_subhalo* first_progenitor_;
    internal_subhalo* next_progenitor_;
    internal_subhalo* descendant_;
    internal_subhalo* first_subhalo_in_fof_group_;
    internal_subhalo* next_subhalo_in_fof_group_;
    internal_subhalo* last_progenitor_;
    internal_subhalo* main_leaf_progenitor_;
    internal_subhalo* root_descendant_;

    /** Constructor. */
    internal_subhalo(const data_format data)
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

  template <typename MAP>
  void create_subhalo_map(const std::string& treefilename, MAP& sub_map) {
    // Read data
    std::cout << "Reading data from file " << treefilename << "\n";
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

    // Create internal_subhalo objects, storing pointers to them in a map.
    std::cout << "Creating subhalo map...\n";
    uint64_t nrows = SubhaloID.size();
    for (uint64_t rownum = 0; rownum < nrows; ++rownum) {
      sub_map.emplace(SubhaloID[rownum], new internal_subhalo(
          data_format(
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
  }

  //////////////////////////////
  // PRIVATE MEMBER VARIABLES //
  //////////////////////////////

  /** (snapnum, subfind_id) to Subhalo* mapping. */
  std::vector<std::vector<internal_subhalo*>> subhalos_;

};
