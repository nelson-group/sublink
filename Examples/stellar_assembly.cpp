/** @file stellar_assembly.cpp
 * @brief For each particle with parttype=4 (mostly stars), determine if
 *        it was formed in-situ or ex-situ, along with other information.
 *
 * @author Vicente Rodriguez-Gomez (vrg@jhu.edu)
 */

// Include some extra quantities from the merger trees:
#define COUNT_MERGERS

// Include main_leaf_progenitor
#define EXTRA_POINTERS

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <iomanip>  // setfill, setw, setprecision
#include <algorithm>  // std::max

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
#include "../Spatial/Point.hpp"

// Type of stellar particles
static constexpr int parttype = 4;

/** @brief Synonym for Tree::Snapshot from ReadTreeHDF5.hpp. */
typedef Tree::Snapshot Snapshot;
/** @brief Synonym for Tree::Subhalo from ReadTreeHDF5.hpp. */
typedef Tree::Subhalo Subhalo;

/** @brief Datatype for stellar particles. */
struct ParticleInfo{
  /** @brief ID of this particle. */
  part_id_type id;
  /** @brief SubfindID of the subhalo where the particle formed. */
  index_type subfind_id_at_formation;
  /** @brief Distance to galactic center at the time of formation */
  real_type distance_at_formation;
  /** @brief SnapNum of the subhalo where the particle formed. */
  snapnum_type snapnum_at_formation;
  /** Constructor. */
  ParticleInfo(part_id_type id_, index_type subf_id_,
      real_type distance_at_formation_, snapnum_type snapnum_)
      : id(id_), subfind_id_at_formation(subf_id_),
        distance_at_formation(distance_at_formation_),
        snapnum_at_formation(snapnum_) {
  }
  /** Disable default constructor. */
  ParticleInfo() = delete;
  /** Custom less-than operator used with std::lower_bound. */
  bool operator<(const part_id_type& other_id) {
    return id < other_id;
  }
};

/** @brief Comparison function to sort by particle ID. */
bool compareByID(const ParticleInfo& a, const ParticleInfo& b) {
  return a.id < b.id;
}

/** @brief Get galactocentric distance of each stellar particle in units of
 * the stellar half-mass radius.
 */
std::vector<real_type> get_galactocentric_distances(
    const std::string& basedir, const snapnum_type snapnum,
    const std::vector<index_type>& SubfindID,
    const std::vector<uint32_t>& sub_offset) {

  uint64_t nparts = SubfindID.size();

  // Get box size from snapshot header
  std::stringstream tmp_stream;
  tmp_stream << basedir << "/snapdir_" <<
      std::setfill('0') << std::setw(3) << snapnum << "/snap_" <<
      std::setfill('0') << std::setw(3) << snapnum << ".0.hdf5";
  std::string file_name = tmp_stream.str();
  auto box_size = static_cast<real_type>(
      arepo::get_scalar_attribute<double>(file_name, "BoxSize"));

  // Read particle positions.
  std::vector<Point> all_pos;
  // Sometimes the coordinates are stored as double type. In that case:
  tmp_stream.str("");
  tmp_stream << "/PartType" << parttype;
  std::string parttype_str = tmp_stream.str();
  auto file = H5::H5File(file_name, H5F_ACC_RDONLY );
  auto group = H5::Group(file.openGroup(parttype_str));
  auto dataset = H5::DataSet(group.openDataSet("Coordinates"));
  auto dt = dataset.getDataType();
  file.close();
  if (dt.getSize() == 8) {
    auto all_pos_double = arepo::read_block<DoubleArray<3>>(
        basedir, snapnum, "Coordinates", parttype);
    std::cout << "NOTE: converting coordinates from double to float...\n";
    all_pos = std::vector<Point>(all_pos_double.begin(),
        all_pos_double.end());
  }
  else { // hopefully dt.getSize() == 4
    all_pos = arepo::read_block<Point>(
        basedir, snapnum, "Coordinates", parttype);
  }
  assert(all_pos.size() == nparts);

  // Read stellar half-mass radius
  auto sub_halfmassrad = subfind::read_block<real_type>(
      basedir, snapnum, "Subhalo", "SubhaloHalfmassRadType", parttype);

  // Store results here
  std::vector<real_type> ParticleDistance(nparts, -1);

  // Iterate over stellar particles.
  for (uint64_t i = 0; i < SubfindID.size(); ++i) {
    auto sub_index = SubfindID[i];
    // Only proceed if stellar particle was actually formed inside a subhalo.
    if (sub_index == -1)
      continue;
    // Take boundary conditions into account
    auto dx = all_pos[i] - all_pos[sub_offset[SubfindID[i]]];
    for (int k = 0; k < 3; ++k) {
      if (std::abs(dx[k]) > 0.5*box_size)
        dx[k] = dx[k] - std::copysign(box_size, dx[k] - 0.5*box_size);
    }
    // Only store if stellar half-mass radius is well defined.
    if (sub_halfmassrad[sub_index] > 0)
      ParticleDistance[i] = norm(dx) / sub_halfmassrad[sub_index];
  }
  return ParticleDistance;
}

/** @brief Populate @a cur_stars with info from given snapshot,
 * which acts as a sort of restart file.
 */
void read_stars(const std::string& writepath,
    const snapnum_type snapnum,
    std::vector<ParticleInfo>& cur_stars) {

  assert(cur_stars.size() == 0);

  // Filename
  std::stringstream tmp_stream;
  tmp_stream << writepath << "_" <<
      std::setfill('0') << std::setw(3) << snapnum << ".hdf5";
  std::string filename = tmp_stream.str();

  // Read data
  auto ParticleID = read_dataset<part_id_type>(filename, "ParticleID");
  auto SubfindIDAtFormation = read_dataset<index_type>(filename, "SubfindIDAtFormation");
  auto DistanceAtFormation = read_dataset<real_type>(filename, "DistanceAtFormation");
  auto SnapNumAtFormation = read_dataset<snapnum_type>(filename, "SnapNumAtFormation");

  // Add to cur_stars structure
  uint64_t nparts = ParticleID.size();
  cur_stars.reserve(cur_stars.size() + nparts);
  for (uint64_t i = 0; i < nparts; ++i) {
    cur_stars.emplace_back(ParticleID[i], SubfindIDAtFormation[i], DistanceAtFormation[i], SnapNumAtFormation[i]);
  }
}

/** @brief Update @a cur_stars structure with information from new snapshot.
 *
 * @param[in,out] cur_stars Vector with stellar particle information.
 * @param[in] ParticleID Particle IDs of particles from new snapshot.
 * @param[in] SubfindID Subfind IDs of particles from new snapshot.
 * @param[in] ParticleDistance Vector with particle distances.
 * @param[in] snapnum Snapshot number of new snapshot.
 */
void update_stars(std::vector<ParticleInfo>& cur_stars,
    const std::vector<part_id_type>& ParticleID,
    const std::vector<index_type>& SubfindID,
    const std::vector<real_type>& ParticleDistance,
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
    stars_aux.emplace_back(ParticleID[i], SubfindID[i], ParticleDistance[i], snapnum);
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

  // Iterate over particles and add info to cur_stars.
  std::cout << "Dealing with particle coincidences...\n";
   wall_clock.start();
   uint64_t disappeared_particles = 0;
   bool skip_one = false;
   for (auto it = stars_aux.begin(); it != stars_aux.end(); ++it) {
     if (skip_one) {
       skip_one = false;
       continue;
     }
     if ((it+1 != stars_aux.end()) && (it->id == (it+1)->id)) {
       // Check that particles have been added in "chronological" order:
       assert(it->snapnum_at_formation <= (it+1)->snapnum_at_formation);
       assert((it+1)->snapnum_at_formation <= snapnum);
       // Now, we handle duplicate ParticleIDs (yes, this happens
       // in IllustrisTNG) in the following way:
       // 1) If the repeated ID is new (snapnum_at_formation == snapnum),
       // we add both particles to cur_stars (so that the ParticleID vector
       // has the right size, etc.), but they are ignored in later snapshots.
       // 2) If the repeated ID is old (snapnum_at_formation < snapnum),
       // we ignore both particles (anyway, this only happens for ten or so
       // particles at very high redshifts).
       if ((it+1)->snapnum_at_formation < snapnum) {
         std::cerr << "WARNING: Ignoring particles with ParticleID " <<
             it->id << " that formed in snapshots " << it->snapnum_at_formation <<
             " and " << (it+1)->snapnum_at_formation << ".\n";
         skip_one = true;
         continue;
       }
       else { // (it+1)->snapnum_at_formation == snapnum
         if (it->snapnum_at_formation == snapnum) {
           std::cerr << "WARNING: ParticleID " << it->id <<
               " appears twice in snapshot " << snapnum << ".\n";
         }
         else { // it->snapnum_at_formation < snapnum
           // All good. Particle is not new. Only add old info to cur_stars.
           skip_one = true;
         }
       }
     }
     else { // it->id != (it+1)->id (or this is the last particle)
       if (it->snapnum_at_formation < snapnum) {
         // Particle disappeared from current snapshot (wind particle?).
         // Do not add to updated cur_stars.
         disappeared_particles += 1;
         continue;
       }
       else {
         // New stellar particle. Add to cur_stars.
         assert(it->snapnum_at_formation == snapnum);
       }
     }
     cur_stars.push_back(*it);
   }
   std::cout << "Time: " << wall_clock.seconds() << " s.\n";
   std::cout << disappeared_particles << " of " << cur_stars.size() <<
       " particles disappeared in this snapshot.\n";
   assert(cur_stars.size() == ParticleID.size());
}

/** @brief Carry out stellar assembly calculations. */
void stellar_assembly(const std::string& basedir, const std::string& treedir,
    const std::string& writepath, const snapnum_type snapnum_first,
    const snapnum_type snapnum_last, const snapnum_type snapnum_restart) {

  // Load merger tree
  WallClock wall_clock;
  std::cout << "Loading merger tree...\n";
  int filenum = -1;  // read from all merger tree files
  std::string name = "tree_extended";  // Full format
  Tree tree(treedir, name, filenum);
  std::cout << "Loaded merger tree. Total time: " <<
      wall_clock.seconds() << " s.\n";
  std::cout << "\n";

  // Store (persistent) stellar particle info in this vector.
  std::vector<ParticleInfo> cur_stars;

  // Iterate over snapshots starting with snapnum_start (!= snapnum_first):
  snapnum_type snapnum_start = snapnum_first;

  // Read restart info if necessary.
  wall_clock.start();
  std::cout << "Loading restart info from snapshot " << snapnum_restart << "...\n";
  if (snapnum_restart >= 0) {
    read_stars(writepath, snapnum_restart, cur_stars);
    // Iterate over snapshots starting with snapnum_restart+1
    snapnum_start = snapnum_restart+1;
  }
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  std::cout << "\n";

  // Iterate over snapshots
  std::cout << "Iterating over snapshots..." << std::endl;
  for (auto snapnum = snapnum_start; snapnum <= snapnum_last; ++snapnum) {

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

    // Calculate galactocentric distances
    std::cout << "Calculating galactocentric distances...\n";
    wall_clock.start();
    auto ParticleDistance = get_galactocentric_distances(basedir,
        snapnum, SubfindID, sub_offset);
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    // Update cur_stars with info from current snapshot.
    update_stars(cur_stars, ParticleID, SubfindID, ParticleDistance, snapnum);

    // Initialize some output arrays
    std::cout << "Initializing some output arrays...\n";
    wall_clock.start();
    std::vector<index_type> SubfindIDAtFormation(nparts, -1);
    std::vector<snapnum_type> SnapNumAtFormation(nparts, -1);
    std::vector<int8_t> InSitu(nparts, -1);
    std::vector<int8_t> AfterInfall(nparts, -1);
    std::vector<int8_t> AccretionOrigin(nparts, -1);
    std::vector<real_type> MergerMassRatio(nparts, -1);
    std::vector<real_type> DistanceAtFormation(nparts, -1);
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
      DistanceAtFormation[pos] = cur_info.distance_at_formation;

      // If current star particle does not currently belong to any subhalo,
      // leave most properties undefined (= -1).
      index_type subfind_id = SubfindID[pos];
      if (subfind_id == -1)
        continue;

      // If current star particle was formed outside of any subhalo,
      // define as ex-situ and leave other properties undefined (= -1).
      if (cur_info.subfind_id_at_formation == -1) {
        InSitu[pos] = 0;
        continue;
      }

      // Determine some properties using the merger trees
      auto cur_sub = tree.subhalo(snapnum, subfind_id);
      auto form_sub = tree.subhalo(cur_info.snapnum_at_formation,
                                   cur_info.subfind_id_at_formation);

      if (along_main_branch(cur_sub, form_sub))
        InSitu[pos] = 1;
      else {  // ex situ
        assert(form_sub.snapnum() < cur_sub.snapnum());
        InSitu[pos] = 0;
        AfterInfall[pos] = static_cast<int>(after_infall(cur_sub, form_sub));
        if (is_descendant(cur_sub, form_sub))
          AccretionOrigin[pos] = 0;  // completed merger
        else if (eventually_merge(cur_sub, form_sub))
          AccretionOrigin[pos] = 1;  // ongoing merger
        else
          AccretionOrigin[pos] = 2;  // "flyby"
        // Merger mass ratio can be determined in any of the above cases:
        MergerMassRatio[pos] = get_merger_mass_ratio(cur_sub, form_sub);
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
    add_array(writefile, AccretionOrigin, "AccretionOrigin",
        H5::PredType::NATIVE_INT8);
    add_array(writefile, MergerMassRatio, "MergerMassRatio",
        H5::PredType::NATIVE_FLOAT);
    add_array(writefile, DistanceAtFormation, "DistanceAtFormation",
        H5::PredType::NATIVE_FLOAT);

    writefile.close();
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    std::cout << "Finished for snapshot " << snapnum << ".\n";
    std::cout << "Total time: " << wall_clock_all.seconds() << " s.\n";
    std::cout << "\n";
  }
}

/** @brief For each particle with parttype=4 (mostly stars), determine if
 *         it was formed in-situ or ex-situ, along with other information.
 */
int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 7) {
    std::cerr << "Usage: " << argv[0] << " basedir treedir writepath" <<
        " snapnum_first snapnum_last snapnum_restart\n";
    exit(1);
  }

  // Read input
  std::string basedir(argv[1]);
  std::string treedir(argv[2]);
  std::string writepath(argv[3]);
  int16_t snapnum_first = atoi(argv[4]);
  int16_t snapnum_last = atoi(argv[5]);
  int16_t snapnum_restart = atoi(argv[6]);

  // Measure CPU and wall clock (real) time
  WallClock wall_clock;
  CPUClock cpu_clock;

  // Do stuff
  stellar_assembly(basedir, treedir, writepath, snapnum_first, snapnum_last, snapnum_restart);

  // Print wall clock time and speedup
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  std::cout << "Speedup: " << cpu_clock.seconds()/wall_clock.seconds() << ".\n";

  return 0;
}
