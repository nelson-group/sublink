/** @file stellar_assembly.cpp
 * @brief Determine the ex situ stellar fraction of galaxies using
 *        different methods.
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <iomanip>  // setfill, setw, setprecision
#include <chrono>   // Wall clock time
#include <ctime>    // CPU time

#include "../InputOutput/ReadArepoHDF5.hpp"
#include "../InputOutput/ReadSubfindHDF5.hpp"
#include "../InputOutput/GeneralHDF5.hpp"
#include "../InputOutput/ReadTreeHDF5.hpp"
#include "../Util/TreeUtil.hpp"

/** @brief Type of particle IDs. */
typedef uint64_t part_id_type;

/** Get some types from ReadTreeHDF5.hpp and make them our own. */
typedef Tree::Snapshot Snapshot;
typedef Tree::Subhalo Subhalo;
typedef typename Tree::sub_id_type sub_id_type;
typedef typename Tree::index_type index_type;
typedef typename Tree::snapnum_type snapnum_type;

// Type containing basic formation info for each stellar particle.
struct formation_info {
  index_type subfind_id_at_formation;
  snapnum_type snapnum_at_formation;

  /** Default constructor. Creates invalid object. */
  formation_info()
      : subfind_id_at_formation(-1), snapnum_at_formation(-1) {
  }
  /** Constructor. */
  formation_info(index_type subf_id_, snapnum_type snapnum_)
      : subfind_id_at_formation(subf_id_), snapnum_at_formation(snapnum_) {
  }
};

/** @brief Carry out stellar assembly calculations.
 *
 *
 * TO DO: http://en.cppreference.com/w/cpp/language/parameter_pack
 */
void stellar_assembly(const std::string& basedir, const std::string& treedir,
    const std::string& writepath, const snapnum_type snapnum_first,
    const snapnum_type snapnum_last) {

  const int parttype = 4;  // stars

  // Load the (SnapNum,SubfindID)->SubhaloID matching for all snapshots
  auto start = std::chrono::system_clock::now();
  std::cout << "Creating (SnapNum,SubfindID)->SubhaloID mapping...\n";
  std::vector<std::vector<sub_id_type>> SubfindIDToSubhaloID;
  for (auto snapnum = snapnum_first; snapnum < snapnum_last+1; ++snapnum) {
    std::stringstream tmp_stream;
    tmp_stream << treedir << "/offsets/offsets_" <<
        std::setfill('0') << std::setw(3) << snapnum << ".hdf5";
    std::string filename_offsets = tmp_stream.str();
    SubfindIDToSubhaloID.push_back(read_dataset<sub_id_type>(
        filename_offsets, "SubhaloID"));
  }
  std::cout << "Time: " << std::chrono::duration<double>(
      std::chrono::system_clock::now()-start).count() << " s.\n";
  std::cout << "\n";

  // Store formation info for each particle in this map:
  std::map<part_id_type, formation_info> ParticleIDMap;

  // Iterate over snapshots
  std::cout << "Iterating over snapshots..." << std::endl;
  for (auto snapnum = snapnum_first; snapnum < snapnum_last+1; ++snapnum) {

    // AD HOC: if L75n1820FP, skip snapshots 53 and 55.
    if (basedir == "/n/hernquistfs1/Illustris/Runs/L75n1820FP/output") {
      if ((snapnum == 53) || (snapnum == 55)) {
        std::cout << "WARNING: Skipping snapshot " << snapnum << ".\n";
        continue;
      }
    }

    // Snapshot filename, without the file number
    std::stringstream tmp_stream;
    tmp_stream << basedir << "/snapdir_" <<
        std::setfill('0') << std::setw(3) << snapnum << "/snap_" <<
        std::setfill('0') << std::setw(3) << snapnum;
    std::string snapname = tmp_stream.str();

    // Initialize output arrays
    start = std::chrono::system_clock::now();
    std::cout << "Initializing arrays...\n";
    auto ParticleID = arepo::read_block<part_id_type>(snapname, "ParticleIDs",
        parttype);
    uint64_t nparts = ParticleID.size();
    std::vector<index_type> SubfindID(nparts, -1);
    std::vector<index_type> SubfindIDAtFormation(nparts, -1);
    std::vector<snapnum_type> SnapNumAtFormation(nparts, -1);
    std::vector<int8_t> InSitu(nparts, -1);
    std::vector<int8_t> BeforeInfall(nparts, -1);
    std::vector<int8_t> BeforeR200(nparts, -1);
    std::cout << "Time: " << std::chrono::duration<double>(
        std::chrono::system_clock::now()-start).count() << " s.\n";

    // Output filename
    tmp_stream.str("");
    tmp_stream << writepath << "_" <<
        std::setfill('0') << std::setw(3) << snapnum << ".hdf5";
    std::string writefilename = tmp_stream.str();

    // Only proceed if there are at least some stellar particles
    std::cout << "There are " << nparts << " stellar particles in snapshot " <<
        snapnum << "." << std::endl;
    if (nparts == 0) {
      // Write to file and continue
      H5::H5File* file = new H5::H5File(writefilename, H5F_ACC_TRUNC);
      add_array(file, ParticleID, "ParticleID",
          H5::PredType::NATIVE_UINT64);
      add_array(file, SubfindID, "SubfindID",
          H5::PredType::NATIVE_INT32);
      add_array(file, SubfindIDAtFormation, "SubfindIDAtFormation",
          H5::PredType::NATIVE_INT32);
      add_array(file, SnapNumAtFormation, "SnapNumAtFormation",
          H5::PredType::NATIVE_INT16);
      add_array(file, InSitu, "InSitu",
          H5::PredType::NATIVE_INT8);
      add_array(file, BeforeInfall, "BeforeInfall",
          H5::PredType::NATIVE_INT8);
      add_array(file, BeforeR200, "BeforeR200",
          H5::PredType::NATIVE_INT8);
      file->close();
      delete file;
      std::cout << "Finished for snapshot " << snapnum << ".\n\n";
      continue;
    }

    // Only proceed if there are at least some subhalos
    auto nsubs = subfind::get_scalar_attribute<uint32_t>(basedir, snapnum,
        "Nsubgroups_Total");
    auto sub_len = subfind::read_block<uint32_t>(basedir, snapnum, "Subhalo",
        "SubhaloLenType", parttype);
    assert(nsubs == sub_len.size());
    std::cout << "There are " << nsubs << " subhalos in snapshot " <<
        snapnum << ".\n";
    if (nsubs == 0) {
      // Write to file and continue
      H5::H5File* file = new H5::H5File(writefilename, H5F_ACC_TRUNC);
      add_array(file, ParticleID, "ParticleID",
          H5::PredType::NATIVE_UINT64);
      add_array(file, SubfindID, "SubfindID",
          H5::PredType::NATIVE_INT32);
      add_array(file, SubfindIDAtFormation, "SubfindIDAtFormation",
          H5::PredType::NATIVE_INT32);
      add_array(file, SnapNumAtFormation, "SnapNumAtFormation",
          H5::PredType::NATIVE_INT16);
      add_array(file, InSitu, "InSitu",
          H5::PredType::NATIVE_INT8);
      add_array(file, BeforeInfall, "BeforeInfall",
          H5::PredType::NATIVE_INT8);
      add_array(file, BeforeR200, "BeforeR200",
          H5::PredType::NATIVE_INT8);
      file->close();
      delete file;
      std::cout << "Finished for snapshot " << snapnum << ".\n\n";
      continue;
    }

    // Associate particles with subhalos in this snapshot
    start = std::chrono::system_clock::now();
    std::cout << "Associating particles with subhalos...\n";
    auto sub_offset = calculate_subhalo_offsets(basedir, snapnum, parttype);
    for (uint32_t sub_uindex = 0; sub_uindex < nsubs; ++sub_uindex) {
      part_id_type snap_count = sub_offset[sub_uindex];
      auto cur_sub_len = sub_len[sub_uindex];
      for (uint32_t i = 0; i < cur_sub_len; ++i) {
        SubfindID[snap_count] = sub_uindex;
        ++snap_count;
      }
    }
    std::cout << "Time: " << std::chrono::duration<double>(
        std::chrono::system_clock::now()-start).count() << " s.\n";

    // We will also need some merger tree info
    tmp_stream.str("");
    tmp_stream << treedir << "/offsets/offsets_" <<
        std::setfill('0') << std::setw(3) << snapnum << ".hdf5";
    std::string filename_offsets = tmp_stream.str();
    auto SubhaloID = read_dataset<sub_id_type>(filename_offsets, "SubhaloID");
    auto MainLeafProgenitorID = read_dataset<sub_id_type>(filename_offsets,
        "MainLeafProgenitorID");

    // Iterate over stellar particles to fill arrays
    start = std::chrono::system_clock::now();
    std::cout << "Iterating over stellar particles...\n";
    for (uint64_t pos = 0; pos < nparts; ++pos) {
      formation_info cur_info;
      auto map_it = ParticleIDMap.find(ParticleID[pos]);
      if (map_it == ParticleIDMap.end()) {
        // Add new stellar particle to map
        cur_info = formation_info(SubfindID[pos], snapnum);
        ParticleIDMap[ParticleID[pos]] = cur_info;
      }
      else
        cur_info = map_it->second;

      SubfindIDAtFormation[pos] = cur_info.subfind_id_at_formation;
      SnapNumAtFormation[pos] = cur_info.snapnum_at_formation;

      // We determine the "in situ" property using the merger trees
      index_type subfind_id = SubfindID[pos];
      sub_id_type subhalo_id_at_formation = -1;
      if (cur_info.subfind_id_at_formation != -1)
        subhalo_id_at_formation = SubfindIDToSubhaloID
            [cur_info.snapnum_at_formation][cur_info.subfind_id_at_formation];

      if ((subfind_id == -1) || (subhalo_id_at_formation == -1)) {
        // If the star particle was formed outside of any subhalo,
        // or is not currently found inside any subhalo,
        // we consider it to be "ex situ."
        InSitu[pos] = 0;
      }
      else {
        if ((subhalo_id_at_formation >= SubhaloID[subfind_id]) &&
            (subhalo_id_at_formation <= MainLeafProgenitorID[subfind_id]))
          InSitu[pos] = 1;
        else {
          InSitu[pos] = 0;

          // Also check if subhalo_id_at_formation and the subhalo
          // from the same snapshot along the main branch are inside
          // the same FoF group

          // Also check if the stellar particle is within R200
          // of the parent FoF group.
        }
      }
    }
    std::cout << "Time: " << std::chrono::duration<double>(
        std::chrono::system_clock::now()-start).count() << " s.\n";

    // NEW: In order to clear some memory, create new map only with
    // particles that currently exist...
    start = std::chrono::system_clock::now();
    std::cout << "Refreshing map...\n";
    std::map<part_id_type, formation_info> ParticleIDMap_aux;
    for (uint64_t pos = 0; pos < nparts; pos++) {
      auto map_it = ParticleIDMap.find(ParticleID[pos]);
      if (map_it == ParticleIDMap.end())
        std::cerr << "Something went wrong: particle should be in map.\n";
      else
        ParticleIDMap_aux[ParticleID[pos]] = map_it->second;
    }
    ParticleIDMap_aux.swap(ParticleIDMap);
    ParticleIDMap_aux.clear();
    std::cout << "Time: " << std::chrono::duration<double>(
        std::chrono::system_clock::now()-start).count() << " s.\n";

    // Write to file.
    start = std::chrono::system_clock::now();
    std::cout << "Writing to file...\n";
    H5::H5File* file = new H5::H5File(writefilename, H5F_ACC_TRUNC);
    add_array(file, ParticleID, "ParticleID",
        H5::PredType::NATIVE_UINT64);
    add_array(file, SubfindID, "SubfindID",
        H5::PredType::NATIVE_INT32);
    add_array(file, SubfindIDAtFormation, "SubfindIDAtFormation",
        H5::PredType::NATIVE_INT32);
    add_array(file, SnapNumAtFormation, "SnapNumAtFormation",
        H5::PredType::NATIVE_INT16);
    add_array(file, InSitu, "InSitu",
        H5::PredType::NATIVE_INT8);
    add_array(file, BeforeInfall, "BeforeInfall",
        H5::PredType::NATIVE_INT8);
    add_array(file, BeforeR200, "BeforeR200",
        H5::PredType::NATIVE_INT8);
    file->close();
    delete file;
    std::cout << "Time: " << std::chrono::duration<double>(
        std::chrono::system_clock::now()-start).count() << " s.\n";

    std::cout << "Finished for snapshot " << snapnum << ".\n\n";
  }
}

int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 6) {
    std::cerr << "Usage: " << argv[0] << " basedir treedir writepath" <<
        " snapnum_first snapnum_last\n";
    exit(1);
  }

  // Read input
  std::string basedir(argv[1]);
  std::string treedir(argv[2]);
  std::string writepath(argv[3]);
  int16_t snapnum_first = atoi(argv[4]);
  int16_t snapnum_last = atoi(argv[5]);

  // Measure CPU and wall clock (real) time
  auto c_start  = std::clock();
  auto t_start = std::chrono::system_clock::now();

  // Do stuff
  stellar_assembly(basedir, treedir, writepath, snapnum_first, snapnum_last);

  // Print CPU and wall clock time
  std::cout << "Finished.\n";
  std::cout << "CPU time: "  <<
      1.0 * (std::clock()-c_start) / CLOCKS_PER_SEC << " s.\n";
  std::cout << "Wall clock time: " << std::chrono::duration<double>(
      std::chrono::system_clock::now()-t_start).count() << " s.\n";

  return 0;
}
