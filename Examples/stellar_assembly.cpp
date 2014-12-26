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
#include <iomanip>  // setfill, setwidth
#include <chrono>

#include "../InputOutput/ReadArepoHDF5.hpp"

/** @brief Type of particle IDs. */
typedef uint64_t part_id_type;
/** @brief Type of subhalo IDs in the merger trees. */
typedef int64_t sub_id_type;
/** @brief Type of subhalo indices in the Subfind catalogs. */
typedef int32_t subf_id_type;
/** @brief Type of snapshot numbers. */
typedef int16_t snapnum_type;
/** @brief Type of indices and sizes. */
typedef uint64_t size_type;

/** @brief Type associated with each particle inside map structure. */
struct particle_info {
  subf_id_type subfind_id_at_formation;
  snapnum_type snapnum_at_formation;
};

/** @brief Carry out stellar assembly calculations.
 *
 */
void stellar_assembly(const std::string& basedir, const std::string& treedir,
    const std::string& writepath, const snapnum_type snapnum_first,
    const snapnum_type snapnum_last) {

  const int parttype = 4;  // stars

//  // Store stellar particle info in these vectors:
//  std::vector<part_id_type> ParticleID;
//  std::vector<subf_id_type> SubfindID;
//  std::vector<subf_id_type> SubfindIDAtFormation;
//  std::vector<snapnum_type> SnapNumAtFormation;
//  std::vector<int8_t> InSitu;
//  std::vector<int8_t> BeforeInfall;
//  std::vector<int8_t> BeforeR200;

  // Store formation info for each particle in this map:
  std::map<part_id_type, particle_info> ParticleIDMap;

  // Other helpful arrays and variables
  std::vector<uint32_t> sub_len, sub_offset;
  uint32_t nsubs, sub_uindex, snap_count, i;
  int32_t subfind_id;

  // Merger tree info
  std::string filename_offsets;
  sub_id_type subhalo_id_at_formation;
  std::vector<sub_id_type> SubhaloID;
  std::vector<sub_id_type> MainLeafProgenitorID;

  // Load the SubfindID->SubhaloID matching for all snapshots
  std::cout << "Loading SubfindID->SubhaloID matching for all snapshots...\n";
  std::vector<std::vector<sub_id_type>> SubfindIDToSubhaloID;
  for (auto snapnum = snapnum_first; snapnum < snapnum_last+1; ++snapnum) {
    std::stringstream tmp_stream;
    tmp_stream << treedir << "/offsets/offsets_" <<
        std::setfill('0') << std::setw(3) << snapnum << ".hdf5";
    filename_offsets = tmp_stream.str();
    SubfindIDToSubhaloID.push_back(read_dataset<sub_id_type>(
        filename_offsets, "SubhaloID"));
  }

  // Iterate over snapshots
  std::cout << "Iterating over snapshots..." << std::endl;
  for (auto snapnum = snapnum_first; snapnum < snapnum_last+1; snapnum++) {

    // AD HOC: if L75n1820FP, skip snapshots 53 and 55.
    if (basedir == "/n/hernquistfs1/Illustris/Runs/L75n1820FP/output") {
      if ((snapnum == 53) || (snapnum == 55)) {
        std::cout << "WARNING: Skipping snapshot " << snapnum << ".\n";
        continue;
      }
    }

    // Output filename
    std::stringstream tmp_stream;
    tmp_stream << writepath << "_" <<
        std::setfill('0') << std::setw(3) << snapnum << ".hdf5";
    std::string writefilename = tmp_stream.str();

    // Initialize the arrays
    std::cout << "Initializing arrays..." << std::endl;
    std::stringstream tmp_stream;
    tmp_stream << basedir << "/snapdir_" <<
        std::setfill('0') << std::setw(3) << snapnum << "/snap_" <<
        std::setfill('0') << std::setw(3) << snapnum;
    std::string snapname = tmp_stream.str();
    auto ParticleID = arepo::read_block<part_id_type>(snapname, "ParticleIDs",
        parttype);
    size_type nparts = ParticleID.size();

    std::vector<subf_id_type> SubfindID(nparts, -1);
    std::vector<subf_id_type> SubfindIDAtFormation(nparts, -1);
    std::vector<snapnum_type> SnapNumAtFormation(nparts, -1);
    std::vector<int8_t> InSitu(nparts, -1);
    std::vector<int8_t> BeforeInfall(nparts, -1);
    std::vector<int8_t> BeforeR200(nparts, -1);

    // Only proceed if there are some stellar particles
    std::cout << "There are " << nparts << " stellar particles in snapshot " <<
        snapnum << "." << std::endl;
    if (nparts == 0) {
      std::cout << std::endl;
      write_to_file(writepath, snapnum, ParticleID, SubfindID,
          SubfindIDAtFormation, SnapNumAtFormation, InSitu);
      continue;
    }

    // Associate particles with subhalos in this snapshot
    std::cout << "Associating particles with subhalos..." << std::endl;
    sub_len = readsubf_read_block<uint32_t>(
        basedir, snapnum, "Subhalo", "SubhaloLenType", parttype);
    nsubs = sub_len.size();

    // Only proceed if there are some subhalos
    std::cout << "There are " << nsubs << " subhalos in snapshot " <<
        snapnum << "." << std::endl;
    if (nsubs == 0) {
      std::cout << std::endl;
      write_to_file(writepath, snapnum, ParticleID, SubfindID,
          SubfindIDAtFormation, SnapNumAtFormation, InSitu);
      continue;
    }

    sub_offset = calculate_subhalo_offsets(basedir, snapnum, parttype);
    snap_count = 0;
    for (sub_uindex = 0; sub_uindex < nsubs; sub_uindex++) {
      snap_count = sub_offset[sub_uindex];
      for (i = 0; i < sub_len[sub_uindex]; i++) {
        SubfindID[snap_count] = sub_uindex;
        snap_count++;
      }
    }

    // We will also need some merger tree info
    tmp_stream.str("");
    tmp_stream << treedir << "/offsets/offsets_" <<
        std::setfill('0') << std::setw(3) << snapnum << ".hdf5";
    filename_offsets = tmp_stream.str();
    SubhaloID = read_block_single_file<int64_t>(
        filename_offsets, "SubhaloID");
    MainLeafProgenitorID = read_block_single_file<int64_t>(
        filename_offsets, "MainLeafProgenitorID");

    // Iterate over star particles to populate arrays
    std::cout << "Iterating over star particles..." << std::endl;
    for (uint64_t pos = 0; pos < nparts; pos++) {
      particle_info cur_info;

      auto map_it = ParticleIDMap.find(ParticleID[pos]);
      if (map_it == ParticleIDMap.end()) {
        // Add new star particle to dictionary
        cur_info.subfind_id_at_formation = SubfindID[pos];
        cur_info.snapnum_at_formation = snapnum;
        ParticleIDMap[ParticleID[pos]] = cur_info;
      }
      else
        cur_info = map_it->second;

      SubfindIDAtFormation[pos] = cur_info.subfind_id_at_formation;
      SnapNumAtFormation[pos] = cur_info.snapnum_at_formation;

      // We determine the "in situ" property using the merger trees
      if (cur_info.subfind_id_at_formation == -1)
        subhalo_id_at_formation = -1;
      else
        subhalo_id_at_formation = SubfindIDToSubhaloID[cur_info.snapnum_at_formation][cur_info.subfind_id_at_formation];
      subfind_id = SubfindID[pos];

      if ((subfind_id == -1) || (subhalo_id_at_formation == -1)) {
        // If the star particle was formed outside of any subhalo,
        // or is not currently found inside any subhalo,
        // we consider this to be "ex situ."
        InSitu[pos] = 0;
      }
      else {
        if ((subhalo_id_at_formation >= SubhaloID[subfind_id]) &&
            (subhalo_id_at_formation <= MainLeafProgenitorID[subfind_id]))
          InSitu[pos] = 1;
        else
          InSitu[pos] = 0;
      }
    }

    // NEW: In order to clear some memory, create new map only with
    // particles that currently exist...
    std::cout << "Refreshing map..." << std::endl;
    std::map<part_id_type, particle_info> ParticleIDMap_aux;
    for (uint64_t pos = 0; pos < nparts; pos++) {
      auto map_it = ParticleIDMap.find(ParticleID[pos]);
      if (map_it == ParticleIDMap.end())
        std::cerr << "Something went wrong: particle should be in map" << std::endl;
      else
        ParticleIDMap_aux[ParticleID[pos]] = map_it->second;
    }

    ParticleIDMap_aux.swap(ParticleIDMap);
    ParticleIDMap_aux.clear();

    // Write to file
    write_to_file(writepath, snapnum, ParticleID, SubfindID,
        SubfindIDAtFormation, SnapNumAtFormation, InSitu);
    std::cout << "Finished for snapshot " << snapnum << "." << std::endl;
  }
}



int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 6) {
    std::std::cerr << "Usage: " << argv[0] << " basedir treedir writepath" <<
        " snapnum_first snapnum_last\n";
    exit(1);
  }

  // Read input
  std::string basedir(argv[1]);
  std::string treedir(argv[2]);
  std::string writepath(argv[3]);
  int16_t snapnum_first = atoi(argv[4]);
  int16_t snapnum_last = atoi(argv[5]);

  // Measure real time
  auto start = std::chrono::system_clock::now();

  // Do stuff
  stellar_assembly(basedir, treedir, writepath, snapnum_first, snapnum_last);

  // Print elapsed time
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Total time: " << elapsed_seconds.count() << " seconds.\n";

  return 0;
}
