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

#include <algorithm>  // stable_sort, lower_bound
#ifdef USE_OPENMP
#include <parallel/algorithm>  // parallel stable_sort
#endif

#include "../InputOutput/ReadArepoHDF5.hpp"
#include "../InputOutput/ReadSubfindHDF5.hpp"
#include "../InputOutput/ReadTreeHDF5.hpp"
#include "../Util/SnapshotUtil.hpp"
#include "../Util/GeneralUtil.hpp"
#include "../Util/TreeUtil.hpp"

// Type of stellar particles
static constexpr int parttype = 4;

/** Get some types from ReadTreeHDF5.hpp and make them our own. */
typedef Tree::Snapshot Snapshot;
typedef Tree::Subhalo Subhalo;

/** Datatype for stellar particles. */
struct ParticleInfo{
  part_id_type id;
  index_type subfind_id_at_formation;
  snapnum_type snapnum_at_formation;
  /** Constructor. */
  ParticleInfo(part_id_type id_, index_type subf_id_, snapnum_type snapnum_)
      : id(id_), subfind_id_at_formation(subf_id_),
        snapnum_at_formation(snapnum_) {
  }
  /** Custom less-than operator used with std::lower_bound. */
  bool operator<(const part_id_type& other_id) {
    return id < other_id;
  }
};

/** Comparison function to sort by particle ID. */
bool compareByID(const ParticleInfo& a, const ParticleInfo& b) {
  return a.id < b.id;
}

/** @brief Update @a cur_stars structure with information from new snapshot.
 *
 * @param[in,out] cur_stars Vector with stellar particle information.
 * @param[in] ParticleID Particle IDs of particles from new snapshot.
 * @param[in] SubfindID Subfind IDs of particles from new snapshot.
 * @param[in] basedir Directory with simulation output.
 * @param[in] snapnum Snapshot number of new snapshot.
 */
void update_stars(std::vector<ParticleInfo>& cur_stars,
    const std::vector<part_id_type>& ParticleID,
    const std::vector<index_type>& SubfindID,
    const snapnum_type snapnum) {

  std::vector<ParticleInfo> stars_aux;
  stars_aux.swap(cur_stars);
  assert(cur_stars.size() == 0);

  // Add stellar particles from new snapshot
  std::cout << "Extending array with particles from new snapshot...\n";
  WallClock wall_clock;
  uint64_t nparts = ParticleID.size();
  stars_aux.reserve(stars_aux.size() + nparts);
  for (uint64_t i = 0; i < nparts; ++i)
    stars_aux.emplace_back(ParticleInfo(ParticleID[i], SubfindID[i], snapnum));
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";

  // Sort array
  std::cout << "Sorting array...\n";
  wall_clock.start();
  CPUClock cpu_clock;
#ifdef USE_OPENMP
  __gnu_parallel::stable_sort(stars_aux.begin(), stars_aux.end(), compareByID);
#else
  std::stable_sort(stars_aux.begin(), stars_aux.end(), compareByID);
#endif
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  std::cout << "Speedup: " << cpu_clock.seconds()/wall_clock.seconds() << ".\n";

  // Iterate over particles and add info to cur_stars, skipping over
  // particles that were formed before the current snapshot.
  std::cout << "Finding particle coincidences...\n";
  wall_clock.start();
  bool skip_one = false;
  for (auto it = stars_aux.begin(); it != stars_aux.end(); ++it) {
    if (skip_one) {
      skip_one = false;
      continue;
    }
    if ((it+1 != stars_aux.end()) && (it->id == (it+1)->id)) {
      // Particle is not new. Keep old info.
      skip_one = true;
    }
    else if (it->snapnum_at_formation != snapnum)
      // Particle disappeared from current snapshot. Skip.
      continue;
    cur_stars.push_back(*it);
  }
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  assert(cur_stars.size() == ParticleID.size());
}

/** @brief Carry out stellar assembly calculations. */
void stellar_assembly(const std::string& basedir, const std::string& treedir,
    const std::string& writepath, const snapnum_type snapnum_first,
    const snapnum_type snapnum_last) {

  // Load merger tree
  WallClock wall_clock;
  std::cout << "Loading merger tree...\n";
  int filenum = -1;  // "concatenated" tree file
  std::string name = "tree_extended";  // Full format
  Tree tree(treedir, name, filenum);
  std::cout << "Loaded merger tree. Total time: " <<
      wall_clock.seconds() << " s.\n";
  std::cout << "\n";

  // Store (persistent) stellar particle info in this vector.
  std::vector<ParticleInfo> cur_stars;

  // Iterate over snapshots
  std::cout << "Iterating over snapshots..." << std::endl;
  for (auto snapnum = snapnum_first; snapnum <= snapnum_last; ++snapnum) {

    // AD HOC: if L75n1820FP, skip snapshots 53 and 55.
    if (basedir == "/n/hernquistfs1/Illustris/Runs/L75n1820FP/output") {
      if ((snapnum == 53) || (snapnum == 55)) {
        std::cout << "WARNING: Skipping snapshot " << snapnum << ".\n";
        continue;
      }
    }

    // For performance checks
    WallClock wall_clock_all;

    // Open output file
    std::stringstream tmp_stream;
    tmp_stream << writepath << "_" <<
        std::setfill('0') << std::setw(3) << snapnum << ".hdf5";
    std::string writefilename = tmp_stream.str();
    H5::H5File writefile(writefilename, H5F_ACC_TRUNC);

    // Only proceed if there is at least one subhalo
    auto nsubs = subfind::get_scalar_attribute<uint32_t>(
        basedir, snapnum, "Nsubgroups_Total");
    if (nsubs == 0) {
      std::cout << "No subhalos in snapshot " << snapnum << ". Skipping...\n";
      writefile.close();
      continue;
    }

    // Read particle IDs
    std::cout << "Reading particle IDs...\n";
    wall_clock.start();
    auto ParticleID = arepo::read_block<part_id_type>(
        basedir, snapnum, "ParticleIDs", parttype);
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    // Only proceed if there is at least one stellar particle
    uint64_t nparts = ParticleID.size();
    if (nparts == 0) {
      std::cout << "No stars in snapshot " << snapnum << ". Skipping...\n";
      writefile.close();
      continue;
    }

    // Get Subfind ID of each stellar particle.
    std::cout << "Getting Subfind IDs of particles...\n";
    wall_clock.start();
    std::vector<index_type> SubfindID(nparts, -1);
    auto sub_len = subfind::read_block<uint32_t>(
        basedir, snapnum, "Subhalo", "SubhaloLenType", parttype);
    auto sub_offset = calculate_subhalo_offsets(basedir, snapnum, parttype);
    for (uint32_t sub_uindex = 0; sub_uindex < nsubs; ++sub_uindex) {
      part_id_type snap_count = sub_offset[sub_uindex];
      auto cur_sub_len = sub_len[sub_uindex];
      for (uint32_t i = 0; i < cur_sub_len; ++i) {
        SubfindID[snap_count] = sub_uindex;
        ++snap_count;
      }
    }
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    // Update cur_stars with info from current snapshot.
    update_stars(cur_stars, ParticleID, SubfindID, snapnum);

    // Initialize some output arrays
    std::cout << "Initializing some output arrays...\n";
    wall_clock.start();
    std::vector<index_type> SubfindIDAtFormation(nparts, -1);
    std::vector<snapnum_type> SnapNumAtFormation(nparts, -1);
    std::vector<int8_t> InSitu(nparts, -1);
    std::vector<int8_t> AfterInfall(nparts, -1);
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    // Iterate over stellar particles
    wall_clock.start();
    std::cout << "Iterating over stellar particles...\n";
    for (uint64_t pos = 0; pos < nparts; ++pos) {
      // Get formation info of current particle.
      auto it = std::lower_bound(cur_stars.begin(), cur_stars.end(),
          ParticleID[pos]);
      assert(it->id == ParticleID[pos]);
      auto cur_info = *it;

      SubfindIDAtFormation[pos] = cur_info.subfind_id_at_formation;
      SnapNumAtFormation[pos] = cur_info.snapnum_at_formation;

      // If current star particle does not currently belong to any subhalo,
      // leave InSitu and AfterInfall properties undefined (= -1).
      index_type subfind_id = SubfindID[pos];
      if (subfind_id == -1)
        continue;

      // If current star particle was formed outside of any subhalo,
      // define as ex situ and leave AfterInfall property undefined (= -1).
      if (cur_info.subfind_id_at_formation == -1) {
        InSitu[pos] = 0;
        continue;
      }

      // Determine InSitu property using the merger trees
      auto cur_sub = tree.subhalo(snapnum, subfind_id);
      auto form_sub = tree.subhalo(cur_info.snapnum_at_formation,
                                   cur_info.subfind_id_at_formation);
      if (along_main_branch(cur_sub, form_sub))
        InSitu[pos] = 1;
      else {
        InSitu[pos] = 0;
        AfterInfall[pos] = static_cast<int>(after_infall(cur_sub, form_sub));
      }
    }
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    // Write to file.
    wall_clock.start();
    std::cout << "Writing to file...\n";
    add_array(writefile, ParticleID, "ParticleID",
        H5::PredType::NATIVE_UINT64);
    add_array(writefile, SubfindID, "SubfindID",
        H5::PredType::NATIVE_INT32);
    add_array(writefile, SubfindIDAtFormation, "SubfindIDAtFormation",
        H5::PredType::NATIVE_INT32);
    add_array(writefile, SnapNumAtFormation, "SnapNumAtFormation",
        H5::PredType::NATIVE_INT16);
    add_array(writefile, InSitu, "InSitu",
        H5::PredType::NATIVE_INT8);
    add_array(writefile, AfterInfall, "AfterInfall",
        H5::PredType::NATIVE_INT8);
    writefile.close();
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    std::cout << "Finished for snapshot " << snapnum << ".\n";
    std::cout << "Total time: " << wall_clock_all.seconds() << " s.\n";
    std::cout << "\n";
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
  WallClock wall_clock;
  CPUClock cpu_clock;

  // Do stuff
  stellar_assembly(basedir, treedir, writepath, snapnum_first, snapnum_last);

  // Print wall clock time and speedup
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  std::cout << "Speedup: " << cpu_clock.seconds()/wall_clock.seconds() << ".\n";

  return 0;
}
