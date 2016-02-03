/** @file count_pairs.cpp
 * @brief Count galaxy close pairs and print their main properties into files.
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

// Include some extra quantities from the merger trees:
#define COUNT_MERGERS
#define EXTRA_POINTERS

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

/** Datatype to store all the pair data. */
struct PairData {
  // General
  std::vector<real_type> proj_sep;  // Projected separation (kpc/h, comoving)
  std::vector<real_type> los_sep;  // Line-of-sight separation (kpc/h, comoving)
  std::vector<real_type> vel_diff;  // Line-of-sight velocity difference (km/s, physical)

  // Properties of descendant (if any; otherwise all -1)
  std::vector<real_type> mstar_0;
  std::vector<real_type> sfr_0;
  std::vector<real_type> delta_0;
  std::vector<index_type> index_0;
  std::vector<snapnum_type> snapnum_merger;

  // "Instantaneous" properties of pair
  std::vector<real_type> mstar_instant_1;
  std::vector<real_type> sfr_instant_1;
  std::vector<real_type> delta_instant_1;
  std::vector<index_type> index_instant_1;
  std::vector<real_type> mstar_instant_2;
  std::vector<real_type> sfr_instant_2;
  std::vector<real_type> delta_instant_2;
  std::vector<index_type> index_instant_2;
  // Skip snapnum_instant because it's given in the filename.

  // Properties at stellar tmax
  std::vector<real_type> mstar_stmax_1;
  std::vector<real_type> sfr_stmax_1;
  std::vector<real_type> delta_stmax_1;
  std::vector<index_type> index_stmax_1;
  std::vector<real_type> mstar_stmax_2;
  std::vector<real_type> sfr_stmax_2;
  std::vector<real_type> delta_stmax_2;
  std::vector<index_type> index_stmax_2;
  std::vector<snapnum_type> snapnum_stmax;

  // Properties at infall
  std::vector<real_type> mstar_infall_1;
  std::vector<real_type> sfr_infall_1;
  std::vector<real_type> delta_infall_1;
  std::vector<index_type> index_infall_1;
  std::vector<real_type> mstar_infall_2;
  std::vector<real_type> sfr_infall_2;
  std::vector<real_type> delta_infall_2;
  std::vector<index_type> index_infall_2;
  std::vector<snapnum_type> snapnum_infall;

  /** Add a new galaxy pair to this structure. */
  void add_pair(
      real_type proj_sep_,
      real_type los_sep_,
      real_type vel_diff_,

      real_type mstar_0_,
      real_type sfr_0_,
      real_type delta_0_,
      index_type index_0_,
      snapnum_type snapnum_merger_,

      real_type mstar_instant_1_,
      real_type sfr_instant_1_,
      real_type delta_instant_1_,
      index_type index_instant_1_,
      real_type mstar_instant_2_,
      real_type sfr_instant_2_,
      real_type delta_instant_2_,
      index_type index_instant_2_,

      real_type mstar_stmax_1_,
      real_type sfr_stmax_1_,
      real_type delta_stmax_1_,
      index_type index_stmax_1_,
      real_type mstar_stmax_2_,
      real_type sfr_stmax_2_,
      real_type delta_stmax_2_,
      index_type index_stmax_2_,
      snapnum_type snapnum_stmax_,

      real_type mstar_infall_1_,
      real_type sfr_infall_1_,
      real_type delta_infall_1_,
      index_type index_infall_1_,
      real_type mstar_infall_2_,
      real_type sfr_infall_2_,
      real_type delta_infall_2_,
      index_type index_infall_2_,
      snapnum_type snapnum_infall_) {

    proj_sep.push_back(proj_sep_);
    los_sep.push_back(los_sep_);
    vel_diff.push_back(vel_diff_);

    mstar_0.push_back(mstar_0_);
    sfr_0.push_back(sfr_0_);
    delta_0.push_back(delta_0_);
    index_0.push_back(index_0_);
    snapnum_merger.push_back(snapnum_merger_);

    mstar_instant_1.push_back(mstar_instant_1_);
    sfr_instant_1.push_back(sfr_instant_1_);
    delta_instant_1.push_back(delta_instant_1_);
    index_instant_1.push_back(index_instant_1_);
    mstar_instant_2.push_back(mstar_instant_2_);
    sfr_instant_2.push_back(sfr_instant_2_);
    delta_instant_2.push_back(delta_instant_2_);
    index_instant_2.push_back(index_instant_2_);

    mstar_stmax_1.push_back(mstar_stmax_1_);
    sfr_stmax_1.push_back(sfr_stmax_1_);
    delta_stmax_1.push_back(delta_stmax_1_);
    index_stmax_1.push_back(index_stmax_1_);
    mstar_stmax_2.push_back(mstar_stmax_2_);
    sfr_stmax_2.push_back(sfr_stmax_2_);
    delta_stmax_2.push_back(delta_stmax_2_);
    index_stmax_2.push_back(index_stmax_2_);
    snapnum_stmax.push_back(snapnum_stmax_);

    mstar_infall_1.push_back(mstar_infall_1_);
    sfr_infall_1.push_back(sfr_infall_1_);
    delta_infall_1.push_back(delta_infall_1_);
    index_infall_1.push_back(index_infall_1_);
    mstar_infall_2.push_back(mstar_infall_2_);
    sfr_infall_2.push_back(sfr_infall_2_);
    delta_infall_2.push_back(delta_infall_2_);
    index_infall_2.push_back(index_infall_2_);
    snapnum_infall.push_back(snapnum_infall_);
  }
  // Automatically generated default constructor should be OK.
};

/** @brief Return projected separation between @a pos1 and @a pos2. */
real_type get_projected_separation(const Point pos1, const Point pos2) {
  auto dr = Point(pos2[0] - pos1[0], pos2[1] - pos1[1], 0);
  return norm(dr);
}

/** @brief Check if the line-of-sight (z-direction) velocity difference
 * between @a center_gal and @a other_gal is below the specified threshold. */
real_type get_velocity_difference(const Point center_pos, const Point center_vel,
    const Point other_pos, const Point other_vel, const real_type cur_H_kpc_h,
    const real_type cur_z) {

  // Peculiar velocity in km/s (no scaling required, according to wiki)
  real_type peculiar_velocity = other_vel[2] - center_vel[2];
  // Relative position along line of sight (z-axis) in kpc/h (note conversion to physical)
  real_type relative_position = (other_pos[2] - center_pos[2]) / (1.0 + cur_z);
  // Hubble flow in km/s (can be negative)
  real_type hubble_flow = relative_position * cur_H_kpc_h;

  return peculiar_velocity + hubble_flow;
}

/** Return a vector with Galaxies within an annulus centered on @a center_pos
 * with radii (in cylindrical coordinates) between rmin and rmax. */
std::vector<Galaxy> annulus_query(const SpaceSearcherType& s,
    const Point center_pos, const real_type rmin_comoving, const real_type rmax_comoving) {

  // Retrieve points inside this bounding box
  Point center_pos_2d = Point(center_pos[0], center_pos[1], 0.0);
  BoundingBox bb(center_pos_2d - Point(rmax_comoving,rmax_comoving,1.0),
                 center_pos_2d + Point(rmax_comoving,rmax_comoving,1.0));

  // Iterate over points in bounding box using SpaceSearcher
  // and keep the ones that are inside the region of interest.
  std::vector<Galaxy> annulus_gals;
  for (auto it = s.begin(bb); it != s.end(bb); ++it) {
    auto cur_gal = *it;
    auto cur_pos = cur_gal.position();
    auto cur_separation = get_projected_separation(center_pos, cur_pos);
    if ((cur_separation > rmin_comoving) && (cur_separation < rmax_comoving)) {
      annulus_gals.push_back(cur_gal);
    }
  }
  return annulus_gals;
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
void count_pairs_sub(const SpaceSearcherType& s,
    const Tree& tree, const snapnum_type snapnum, const Subhalo center_sub,
    const Point center_pos, const Point center_vel,
    const real_type rmin_comoving, const real_type rmax_comoving,
    const real_type velocity_threshold, const real_type mstar_min,
    const real_type cur_H_kpc_h, const real_type cur_z,
    const std::vector<std::vector<real_type>>& overdensities,
    PairData& pd) {

  assert(center_sub.is_valid());

  // Only proceed if current galaxy is above minimum mass
  if (center_sub.data().SubhaloMassType[parttype_stars] < mstar_min)
    return;

  // Iterate over neighbors (of any mass)
  auto neighbors = annulus_query(s, center_pos, rmin_comoving, rmax_comoving);
  for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
    auto other_gal = *it;  // "spatial" object

    // Just in case (I think I can remove this check)
    if (center_sub.index() == other_gal.value().subfind_id) {
      std::cout << "Skipping repeated galaxy...\n";
      continue;
    }

    // ------------------ CHECK IF VALID CLOSE PAIR -------------------

    // Only proceed if (physical) velocity difference is less than threshold *OR*
    // if line-of-sight separation is less than rmax_comoving
    real_type los_sep = other_gal.position()[2] - center_pos[2];
    real_type vel_diff = get_velocity_difference(center_pos, center_vel,
        other_gal.position(), other_gal.value().velocity, cur_H_kpc_h, cur_z);
    if ((std::abs(vel_diff) > velocity_threshold) &&
        (std::abs(los_sep) > rmax_comoving))
      continue;

    // If we got here, we have found a close pair (either by the vel_diff or
    // los_sep criterion).
    auto other_sub = tree.subhalo(snapnum, other_gal.value().subfind_id);
    assert(other_sub.is_valid());

    // Avoid double counting by making sure that mstar_1 > mstar_2
    // (check merger mass ratio later, e.g. when plotting).
    real_type mstar_instant_1 = center_sub.data().SubhaloMassType[parttype_stars];
    real_type mstar_instant_2 = other_sub.data().SubhaloMassType[parttype_stars];
    if (mstar_instant_1 <= mstar_instant_2)
      continue;

    // ----------------- GET PROPERTIES OF CLOSE PAIR -------------------

    // General properties of close pair
    real_type proj_sep = get_projected_separation(center_pos, other_gal.position());
    // already have los_sep
    // already have vel_diff

    // Properties of descendant (if any)
    real_type mstar_0 = -1;
    real_type sfr_0 = -1;
    real_type delta_0 = -1;
    index_type index_0 = -1;
    snapnum_type snapnum_merger = -1;
    auto desc = get_merger_remnant(center_sub, other_sub);
    if (desc.is_valid()) {
      mstar_0 = desc.data().SubhaloMassType[parttype_stars];
      sfr_0 = desc.data().SubhaloSFR;
      delta_0 = overdensities[desc.snapnum()][desc.index()];
      index_0 = desc.index();
      snapnum_merger = desc.snapnum();
    }

    // Pair properties at current time
    // already have mstar_instant_1
    real_type sfr_instant_1 = center_sub.data().SubhaloSFR;
    real_type delta_instant_1 = overdensities[center_sub.snapnum()][center_sub.index()];
    index_type index_instant_1 = center_sub.index();
    // already have mstar_instant_2
    real_type sfr_instant_2 = other_sub.data().SubhaloSFR;
    real_type delta_instant_2 = overdensities[other_sub.snapnum()][other_sub.index()];
    index_type index_instant_2 = other_sub.index();

    // Properties at (stellar) tmax
    real_type mstar_stmax_1 = -1;
    real_type sfr_stmax_1 = -1;
    real_type delta_stmax_1 = -1;
    index_type index_stmax_1 = -1;
    real_type mstar_stmax_2 = -1;
    real_type sfr_stmax_2 = -1;
    real_type delta_stmax_2 = -1;
    index_type index_stmax_2 = -1;
    snapnum_type snapnum_stmax = -1;
    auto stmax_pair = get_stmax_pair(center_sub, other_sub);
    if (stmax_pair.first.is_valid()) {
      mstar_stmax_1 = stmax_pair.first.data().SubhaloMassType[parttype_stars];
      sfr_stmax_1 = stmax_pair.first.data().SubhaloSFR;
      delta_stmax_1 = overdensities[stmax_pair.first.data().SnapNum][stmax_pair.first.index()];
      index_stmax_1 = stmax_pair.first.index();
      mstar_stmax_2 = stmax_pair.second.data().SubhaloMassType[parttype_stars];
      sfr_stmax_2 = stmax_pair.second.data().SubhaloSFR;
      delta_stmax_2 = overdensities[stmax_pair.second.data().SnapNum][stmax_pair.second.index()];
      index_stmax_2 = stmax_pair.second.index();
      snapnum_stmax = stmax_pair.second.data().SnapNum;
    }

    // Properties at infall
    real_type mstar_infall_1 = -1;
    real_type sfr_infall_1 = -1;
    real_type delta_infall_1 = -1;
    index_type index_infall_1 = -1;
    real_type mstar_infall_2 = -1;
    real_type sfr_infall_2 = -1;
    real_type delta_infall_2 = -1;
    index_type index_infall_2 = -1;
    snapnum_type snapnum_infall = -1;
    auto infall_pair = get_infall_pair(center_sub, other_sub);
    if (infall_pair.first.is_valid()) {
      mstar_infall_1 = infall_pair.first.data().SubhaloMassType[parttype_stars];
      sfr_infall_1 = infall_pair.first.data().SubhaloSFR;
      delta_infall_1 = overdensities[infall_pair.first.data().SnapNum][infall_pair.first.index()];
      index_infall_1 = infall_pair.first.index();
      mstar_infall_2 = infall_pair.second.data().SubhaloMassType[parttype_stars];
      sfr_infall_2 = infall_pair.second.data().SubhaloSFR;
      delta_infall_2 = overdensities[infall_pair.second.data().SnapNum][infall_pair.second.data().SubfindID];
      index_infall_2 = infall_pair.second.index();
      snapnum_infall = infall_pair.second.data().SnapNum;
    }

    // Store pair info in data structure
    pd.add_pair(proj_sep, los_sep, vel_diff, mstar_0, sfr_0, delta_0, index_0,
        snapnum_merger, mstar_instant_1, sfr_instant_1, delta_instant_1, index_instant_1,
        mstar_instant_2, sfr_instant_2, delta_instant_2, index_instant_2, mstar_stmax_1,
        sfr_stmax_1, delta_stmax_1, index_stmax_1, mstar_stmax_2, sfr_stmax_2,
        delta_stmax_2, index_stmax_2, snapnum_stmax, mstar_infall_1, sfr_infall_1,
        delta_infall_1, index_infall_1, mstar_infall_2, sfr_infall_2, delta_infall_2,
        index_infall_2, snapnum_infall);
  }
}

/** @brief Count close pairs and print to files. */
void count_pairs_all(const std::string& simdir, const std::string& treedir,
    const std::string& envdir, const std::string& writepath,
    const snapnum_type snapnum_first, const snapnum_type snapnum_last,
    const real_type rmin, const real_type rmax,
    const real_type velocity_threshold, const real_type log_mstar_min) {

  // Convert minimum mass to 10^10 Msun/h:
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
        std::setfill('0') << std::setw(3) << snapnum << ".hdf5";
    std::string writefilename = tmp_stream.str();

    // Open output HDF5 file
    H5::H5File writefile(writefilename, H5F_ACC_TRUNC);

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
      writefile.close();
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
    // (careful: do not try to access galaxies via AllGalaxies::galaxy(i)
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
    // but probably none of this is necessary because of slow reading speeds).
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

    // Create a structure for storing pair data.
    auto pd = PairData();
    assert(pd.proj_sep.size() == 0);  // just to check

    // Iterate over real (non-ghost) galaxies and count pairs
    std::cout << "Iterating over galaxies...\n";
    wall_clock.start();
    for (auto sub_it = snap.begin(); sub_it != snap.end(); ++sub_it) {
      auto center_sub = *sub_it;
      auto center_pos = sub_pos[center_sub.index()];
      auto center_vel = sub_vel[center_sub.index()];
      count_pairs_sub(s, tree, snapnum, center_sub, center_pos, center_vel,
          rmin_comoving, rmax_comoving, velocity_threshold, mstar_min,
          cur_H_kpc_h, cur_z, overdensities, pd);
    }
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    // Write to HDF5 file.
    add_array(writefile, pd.proj_sep, "proj_sep", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.los_sep, "los_sep", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.vel_diff, "vel_diff", H5::PredType::NATIVE_FLOAT);

    add_array(writefile, pd.mstar_0, "mstar_0", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.sfr_0, "sfr_0", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.delta_0, "delta_0", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.index_0, "index_0", H5::PredType::NATIVE_INT32);
    add_array(writefile, pd.snapnum_merger, "snapnum_merger", H5::PredType::NATIVE_INT16);

    add_array(writefile, pd.mstar_instant_1, "mstar_instant_1", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.sfr_instant_1, "sfr_instant_1", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.delta_instant_1, "delta_instant_1", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.index_instant_1, "index_instant_1", H5::PredType::NATIVE_INT32);
    add_array(writefile, pd.mstar_instant_2, "mstar_instant_2", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.sfr_instant_2, "sfr_instant_2", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.delta_instant_2, "delta_instant_2", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.index_instant_2, "index_instant_2", H5::PredType::NATIVE_INT32);

    add_array(writefile, pd.mstar_stmax_1, "mstar_stmax_1", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.sfr_stmax_1, "sfr_stmax_1", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.delta_stmax_1, "delta_stmax_1", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.index_stmax_1, "index_stmax_1", H5::PredType::NATIVE_INT32);
    add_array(writefile, pd.mstar_stmax_2, "mstar_stmax_2", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.sfr_stmax_2, "sfr_stmax_2", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.delta_stmax_2, "delta_stmax_2", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.index_stmax_2, "index_stmax_2", H5::PredType::NATIVE_INT32);
    add_array(writefile, pd.snapnum_stmax, "snapnum_stmax", H5::PredType::NATIVE_INT16);

    add_array(writefile, pd.mstar_infall_1, "mstar_infall_1", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.sfr_infall_1, "sfr_infall_1", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.delta_infall_1, "delta_infall_1", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.index_infall_1, "index_infall_1", H5::PredType::NATIVE_INT32);
    add_array(writefile, pd.mstar_infall_2, "mstar_infall_2", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.sfr_infall_2, "sfr_infall_2", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.delta_infall_2, "delta_infall_2", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, pd.index_infall_2, "index_infall_2", H5::PredType::NATIVE_INT32);
    add_array(writefile, pd.snapnum_infall, "snapnum_infall", H5::PredType::NATIVE_INT16);

    // Close (and flush) file
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
