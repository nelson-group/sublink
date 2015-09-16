/** @file count_pairs.cpp
 * @brief Count galaxy close pairs and print their main properties into files.
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

// Include some extra quantities from the merger trees:
#define COUNT_MERGERS

#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "../InputOutput/ReadArepoHDF5.hpp"
#include "../InputOutput/ReadSubfindHDF5.hpp"
#include "../InputOutput/ReadTreeHDF5.hpp"
#include "../InputOutput/GeneralHDF5.hpp"
#include "../Util/GeneralUtil.hpp"
#include "../Util/TreeUtil.hpp"
#include "../Util/Cosmology.hpp"  // get_redshifts(), hubble parameter
#include "../Spatial/SpaceSearcher.hpp"
#include "../Spatial/AllGalaxies.hpp"

// Constants
static constexpr int parttype_stars = 4;

/** Custom data structure for galaxy values. */
struct GalaxyData {
  index_type subfind_id;
  Point velocity;

  // No default constructor.
  GalaxyData() = delete;
  /** Constructor. */
  GalaxyData(const index_type subfind_id_, const Point velocity_)
      : subfind_id(subfind_id_), velocity(velocity_) {
  }
};

/** Simple mapping from a Galaxy to its corresponding point on the xy-plane,
 * i.e., making z = 0. */
struct GalaxyToPoint {
  template <typename GAL>
  Point operator()(const GAL& g) {
    auto pos = g.position();
    return Point(pos[0], pos[1], 0.0);
  }
};

// Define Galaxy and AllGalaxies types
// (CAREFUL: do not confuse merger tree stuff, e.g. Subhalo, with
// spatial stuff, e.g. Galaxy)
typedef AllGalaxies<GalaxyData> AllGalaxiesType;
typedef typename AllGalaxiesType::galaxy_type Galaxy;
// Define spatial searcher type
typedef SpaceSearcher<Galaxy, GalaxyToPoint> SpaceSearcherType;

/** @brief Return projected separation between @a pos1 and @a pos2. */
real_type get_projected_separation(const Point pos1, const Point pos2) {
  auto dr = Point(pos2[0] - pos1[0], pos2[1] - pos1[1], 0);
  return norm(dr);
}

/** Return a vector with Galaxies within an annulus centered on @a center_gal
 * with radii (in cylindrical coordinates) rmin and rmax. */
std::vector<Galaxy> annulus_query(const SpaceSearcherType& s,
    const Galaxy center_gal, const real_type rmin, const real_type rmax,
    const real_type box_size) {

  auto center_pos = center_gal.position();

  // Find points inside this bounding box
  BoundingBox bb(center_gal.position() - Point(rmax,rmax,box_size),
                 center_gal.position() + Point(rmax,rmax,box_size));

  // Iterate over points in bounding box using SpaceSearcher
  // and keep the ones that are inside the region of interest.
  std::vector<Galaxy> annulus_gals;
  for (auto it = s.begin(bb); it != s.end(bb); ++it) {
    auto cur_gal = *it;
    auto cur_pos = cur_gal.position();
    auto cur_separation = get_projected_separation(center_pos, cur_pos);
    if ((cur_separation > rmin) && (cur_separation < rmax)) {
      annulus_gals.push_back(cur_gal);
    }
  }
  return annulus_gals;
}

/** @brief Check if the line-of-sight (z-direction) velocity difference
 * between @a center_gal and @a other_gal is below the specified threshold. */
real_type get_velocity_difference(const Galaxy center_gal, const Galaxy other_gal,
    const real_type cur_H_kpc_h) {

  // Relative position along line of sight (z-axis) in kpc/h
  real_type relative_position = other_gal.position()[2] - center_gal.position()[2];
  // Peculiar velocity in km/s (divide by scale factor "a" or something?)
  real_type peculiar_velocity = (other_gal.value().velocity[2] - center_gal.value().velocity[2]);
  // Hubble flow velocity in km/s
  real_type hubble_flow = relative_position * cur_H_kpc_h;  // can be negative

  return peculiar_velocity + hubble_flow;
}

/** @brief Return a vector of vectors with overdensities for all snapshots. */
std::vector<std::vector<real_type>> read_overdensities(
    const std::string& simdir, const std::string& envdir,
    const snapnum_type snapnum_first, const snapnum_type snapnum_last) {

  std::vector<std::vector<real_type>> overdensities(snapnum_last+1);
  for (auto snapnum = snapnum_first; snapnum <= snapnum_last; ++snapnum) {
    // Only proceed if there is at least one subhalo
    auto nsubs = subfind::get_scalar_attribute<uint32_t>(
        simdir + "/output", snapnum, "Nsubgroups_Total");
    if (nsubs == 0)
      continue;

    // AD HOC: if L75n1820FP, skip snapshots 53 and 55.
    if (simdir == "/n/hernquistfs1/Illustris/Runs/L75n1820FP") {
      if ((snapnum == 53) || (snapnum == 55)) {
        std::cout << "WARNING: Skipping snapshot " << snapnum << ".\n";
        continue;
      }
    }

    // Create filename
    std::stringstream tmp_stream;
    tmp_stream << envdir << "/environment_" <<
        std::setfill('0') << std::setw(3) << snapnum << ".hdf5";
    std::string filename = tmp_stream.str();

    // If environment file is empty (it might be empty if there are
    // no "bright" galaxies with r-band magnitude < -19.5),
    // fill array with -999's.
    H5::H5File tmp_file(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (!H5Lexists(tmp_file.getId(), "/delta", H5P_DEFAULT)) {
      std::cout << "Warning: no overdensities for snapshot " << snapnum << ".\n";
      overdensities[snapnum].resize(nsubs, -999);
      tmp_file.close();
      continue;
    }
    tmp_file.close();

    // Read overdensities
    overdensities[snapnum] = read_dataset<real_type>(filename, "delta");
    assert(overdensities[snapnum].size() == nsubs);
  }
  return overdensities;
}

/** @brief Count close pairs for a given subhalo and print data to file. */
void count_pairs_sub(const AllGalaxiesType& ag, const SpaceSearcherType& s,
    const Tree& tree, const snapnum_type snapnum, const Subhalo center_sub,
    const real_type rmin_comoving, const real_type rmax_comoving,
    const real_type velocity_threshold, const real_type mstar_min,
    const real_type cur_H_kpc_h, const real_type box_size,
    const std::vector<std::vector<real_type>>& overdensities,
    std::ofstream& writefile) {

  assert(center_sub.is_valid());

  auto center_gal = ag.galaxy(center_sub.index());  // "spatial" object

  // Only proceed if current galaxy is above minimum mass
  if (center_sub.data().SubhaloMassType[parttype_stars] < mstar_min)
    return;

  // Iterate over neighbors
  auto neighbors = annulus_query(s, center_gal, rmin_comoving, rmax_comoving, box_size);
  for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
    auto other_gal = *it;

    // Just in case
    if (center_gal == other_gal) {
      std::cout << "Skipping repeated galaxy...\n";
      continue;
    }

    // Only proceed if (physical) velocity difference is less than threshold
    real_type vel_diff = get_velocity_difference(center_gal, other_gal, cur_H_kpc_h);
    if (std::abs(vel_diff) > velocity_threshold)
      continue;

    // If we got here, we have found a close pair. Get masses, etc.
    auto other_sub = tree.subhalo(snapnum, other_gal.value().subfind_id);

    assert(other_sub.is_valid());

    // Properties at current snapshot
    writefile << std::setprecision(10) <<
        get_projected_separation(center_gal.position(), other_gal.position()) << "," <<
        other_gal.position()[2] - center_gal.position()[2] << "," <<
        vel_diff << "," <<
        center_sub.data().SubhaloMassType[parttype_stars] << "," <<
        center_sub.data().SubhaloSFR << "," <<
        overdensities[center_sub.snapnum()][center_sub.index()] << "," <<
        other_sub.data().SubhaloMassType[parttype_stars] << "," <<
        other_sub.data().SubhaloSFR << "," <<
        overdensities[other_sub.snapnum()][other_sub.index()] << ",";

    // Properties at (stellar) tmax
    real_type mstar_stmax_1 = -1;
    real_type sfr_stmax_1 = -1;
    real_type overdensity_stmax_1 = -1;
    real_type mstar_stmax_2 = -1;
    real_type sfr_stmax_2 = -1;
    real_type overdensity_stmax_2 = -1;
    snapnum_type snapnum_stmax = -1;
    auto stmax_pair = get_stmax_pair(center_sub, other_sub);
    if (stmax_pair.first.is_valid()) {
      mstar_stmax_1 = stmax_pair.first.data().SubhaloMassType[parttype_stars];
      sfr_stmax_1 = stmax_pair.first.data().SubhaloSFR;
      overdensity_stmax_1 = overdensities[stmax_pair.first.data().SnapNum][stmax_pair.first.data().SubfindID];
      mstar_stmax_2 = stmax_pair.second.data().SubhaloMassType[parttype_stars];
      sfr_stmax_2 = stmax_pair.second.data().SubhaloSFR;
      overdensity_stmax_2 = overdensities[stmax_pair.second.data().SnapNum][stmax_pair.second.data().SubfindID];
      snapnum_stmax = stmax_pair.second.data().SnapNum;
    }
    writefile << std::setprecision(10) <<
        mstar_stmax_1 << "," <<
        sfr_stmax_1 << "," <<
        overdensity_stmax_1 << "," <<
        mstar_stmax_2 << "," <<
        sfr_stmax_2 << "," <<
        overdensity_stmax_2 << "," <<
        snapnum_stmax << ",";

    // Properties at infall
    real_type mstar_infall_1 = -1;
    real_type sfr_infall_1 = -1;
    real_type overdensity_infall_1 = -1;
    real_type mstar_infall_2 = -1;
    real_type sfr_infall_2 = -1;
    real_type overdensity_infall_2 = -1;
    snapnum_type snapnum_infall = -1;
    auto infall_pair = get_infall_pair(center_sub, other_sub);
    if (infall_pair.first.is_valid()) {
      mstar_infall_1 = infall_pair.first.data().SubhaloMassType[parttype_stars];
      sfr_infall_1 = infall_pair.first.data().SubhaloSFR;
      overdensity_infall_1 = overdensities[infall_pair.first.data().SnapNum][infall_pair.first.data().SubfindID];
      mstar_infall_2 = infall_pair.second.data().SubhaloMassType[parttype_stars];
      sfr_infall_2 = infall_pair.second.data().SubhaloSFR;
      overdensity_infall_2 = overdensities[infall_pair.second.data().SnapNum][infall_pair.second.data().SubfindID];
      snapnum_infall = infall_pair.second.data().SnapNum;
    }
    writefile << std::setprecision(10) <<
        mstar_infall_1 << "," <<
        sfr_infall_1 << "," <<
        overdensity_infall_1 << "," <<
        mstar_infall_2 << "," <<
        sfr_infall_2 << "," <<
        overdensity_infall_2 << "," <<
        snapnum_infall << "\n";
  }
}


/** @brief Count close pairs and print to files. */
void count_pairs_all(const std::string& simdir, const std::string& treedir,
    const std::string& envdir, const std::string& writepath,
    const snapnum_type snapnum_first, const snapnum_type snapnum_last,
    const real_type rmin, const real_type rmax,
    const real_type velocity_threshold, const real_type log_mstar_min) {

  // Convert mass to 10^10 Msun/h:
  float mstar_min = cosmo::h * std::pow(10.0, log_mstar_min - 10.0);

  // Get box size in ckpc/h; note conversion to float
  std::stringstream tmp_stream;
  tmp_stream << simdir << "/output/snapdir_" <<
      std::setfill('0') << std::setw(3) << snapnum_last << "/snap_" <<
      std::setfill('0') << std::setw(3) << snapnum_last << ".0.hdf5";
  std::string file_name = tmp_stream.str();
  real_type box_size = static_cast<real_type>(
      arepo::get_scalar_attribute<double>(file_name, "BoxSize"));

  // Get snapshot redshifts.
  auto redshifts_all = cosmo::get_redshifts();

  // Read overdensities for all snapshots
  std::cout << "Reading overdensities...\n";
  WallClock wall_clock;
  auto overdensities = read_overdensities(simdir, envdir,
      snapnum_first, snapnum_last);
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";

  // Load merger tree
  wall_clock.start();
  std::cout << "Loading merger tree...\n";
  int filenum = -1;  // "concatenated" tree file
  std::string name = "tree_extended";  // Full format
  Tree tree(treedir, name, filenum);
  std::cout << "Loaded merger tree. Total time: " <<
      wall_clock.seconds() << " s.\n";
  std::cout << "\n";

  // Iterate over snapshots
  std::cout << "Iterating over snapshots..." << std::endl;
  for (auto snapnum = snapnum_first; snapnum <= snapnum_last; ++snapnum) {
    auto snap = tree.snapshot(snapnum);

    // Some parameters
    real_type cur_z = redshifts_all[snapnum];
    real_type cur_H_kpc_h = cosmo::H_kpc_h(cur_z);
    real_type rmin_comoving = rmin * (1.0 + cur_z);
    real_type rmax_comoving = rmax * (1.0 + cur_z);

    // Filename for this snapshot.
    std::stringstream tmp_stream;
    tmp_stream << writepath << "_" <<
        std::setfill('0') << std::setw(3) << snapnum;
    std::string writefilename = tmp_stream.str();

    // Open file and make sure it's open.
    std::ofstream writefile(writefilename.data());
    if (!writefile.is_open()) {
      std::cerr << "Error: Unable to open file " << writefilename << '\n';
      continue;
    }

    // Load galaxy positions and stellar masses.
    std::cout << "Loading galaxy data...\n";
    WallClock wall_clock;
    std::string basedir = simdir + "/output";
    auto sub_pos = subfind::read_block<Point>(
        basedir, snapnum, "Subhalo", "SubhaloPos", -1);
    auto sub_vel = subfind::read_block<Point>(
        basedir, snapnum, "Subhalo", "SubhaloVel", -1);
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    // Only proceed if nsubs > 0
    uint32_t nsubs = sub_pos.size();
    if (nsubs == 0) {
      std::cout << "Skipping empty snapshot: " << snapnum << ".\n";
      continue;
    }

    // AD HOC: if L75n1820FP, skip snapshots 53 and 55.
    if (simdir == "/n/hernquistfs1/Illustris/Runs/L75n1820FP") {
      if ((snapnum == 53) || (snapnum == 55)) {
        std::cout << "WARNING: Skipping snapshot " << snapnum << ".\n";
        continue;
      }
    }

    // Create Galaxy objects and add to AllGalaxies container
    // (only consider objects that exist in the merger trees).
    std::cout << "Creating galaxy objects...\n";
    wall_clock.start();
    AllGalaxiesType ag;
    for (auto sub_it = snap.begin(); sub_it != snap.end(); ++sub_it) {
      auto sub_index = (*sub_it).index();
      ag.add_galaxy(sub_pos[sub_index], GalaxyData(sub_index, sub_vel[sub_index]));
    }
    // Add one "ghost" galaxy for each dimension
    // (careful: do not try to access ghost galaxies via AllGalaxies::galaxy(i)
    // (optimization: only include ghost galaxies within 100 kpc/h of each border)
    // (Another possible optimization: implement 2D spatial search,
    // but probably none of this is necessary because of slow reading speed).
    for (auto sub_it = snap.begin(); sub_it != snap.end(); ++sub_it) {
      auto sub_index = (*sub_it).index();
      for (int i = 0; i < 3; ++i) {
        Point ghost_pos = sub_pos[sub_index];
        ghost_pos[i] = sub_pos[sub_index][i] - std::copysign(
            box_size, sub_pos[sub_index][i] - 0.5*box_size);
        ag.add_galaxy(ghost_pos, GalaxyData(sub_index, sub_vel[sub_index]));
      }
    }
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    // Create spatial searcher instance.
    std::cout << "Creating SpaceSearcher instance...\n";
    wall_clock.start();
    SpaceSearcherType s(ag.begin(), ag.end(), GalaxyToPoint());
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    // Iterate over real (non-ghost) galaxies and count pairs
    std::cout << "Iterating over galaxies...\n";
    wall_clock.start();
    for (auto sub_it = snap.begin(); sub_it != snap.end(); ++sub_it) {
      auto center_sub = *sub_it;
      count_pairs_sub(ag, s, tree, snapnum, center_sub,
          rmin_comoving, rmax_comoving, velocity_threshold, mstar_min,
          cur_H_kpc_h, box_size, overdensities, writefile);
    }
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    // Flush and close file
    writefile.flush();
    writefile.close();
    std::cout << "Finished for snapshot " << snapnum << ".\n\n";
  }
}


int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 11) {
    std::cerr << "Usage: " << argv[0] << " simdir treedir envdir writepath" <<
        " snapnum_first snapnum_last rmin rmax velocity_threshold log_mstar_min\n";
    exit(1);
  }

  // Read input
  std::string simdir(argv[1]);
  std::string treedir(argv[2]);
  std::string envdir(argv[3]);
  std::string writepath(argv[4]);
  snapnum_type snapnum_first = atoi(argv[5]);
  snapnum_type snapnum_last = atoi(argv[6]);
  real_type rmin = atof(argv[7]);  // ckpc/h
  real_type rmax = atof(argv[8]);  // ckpc/h
  real_type velocity_threshold = atof(argv[9]);  // km/s
  real_type log_mstar_min = atof(argv[10]);  // base 10, Msun

  // Measure CPU and wall clock (real) time
  WallClock wall_clock;
  CPUClock cpu_clock;

  // Do stuff
  count_pairs_all(simdir, treedir, envdir, writepath, snapnum_first, snapnum_last,
      rmin, rmax, velocity_threshold, log_mstar_min);

  // Print wall clock time and speedup
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  std::cout << "Speedup: " << cpu_clock.seconds()/wall_clock.seconds() << ".\n";

  return 0;
}
