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

/** Return a vector with Galaxies within an annulus centered on @a center_gal
 * with radii (in cylindrical coordinates) rmin and rmax. */
std::vector<Galaxy> annulus_query(const SpaceSearcherType& s,
    Galaxy center_gal, real_type rmin, real_type rmax) {

  auto center_pos = center_gal.position();

  // Find points inside this bounding box
  BoundingBox bb(center_gal.position() - Point(rmax,rmax,1),
                 center_gal.position() + Point(rmax,rmax,1));

  // Iterate over points in bounding box using SpaceSearcher
  // and keep the ones that are inside the region of interest.
  std::vector<Galaxy> annulus_gals;
  for (auto it = s.begin(bb); it != s.end(bb); ++it) {
    auto cur_gal = *it;
    auto cur_pos = cur_gal.position();

    auto dr = Point(cur_pos[0] - center_pos[0], cur_pos[1] - center_pos[1], 0);
    if ((norm(dr) > rmin) && (norm(dr) < rmax)) {
      annulus_gals.push_back(cur_gal);
    }
  }

  return annulus_gals;
}


bool within_velocity_threshold(const Galaxy center_gal, const Galaxy other_gal,
    const real_type velocity_threshold, const real_type z) {

  // Relative position along line of sight (z-axis) in kpc/h
  real_type relative_position = other_gal.position()[2] - center_gal.position()[2];
  // Peculiar velocity in km/s (divide by scale factor "a" or something?)
  real_type peculiar_velocity = (other_gal.value().velocity[2] -
      center_gal.value().velocity[2]);
  // Hubble flow velocity in km/s
  real_type hubble_flow = relative_position * cosmo::H_kpc_h(z);  // can be negative

  std::cout << "Relative position: " << relative_position << "kpc/h\n";
  std::cout << "Peculiar velocity: " << peculiar_velocity << "km/s\n";
  std::cout << "Hubble flow: " << hubble_flow << "km/s\n";
  std::cout << "\n";

  return std::abs(peculiar_velocity + hubble_flow) < velocity_threshold;
}


/** @brief Count close pairs and print to files. */
void count_pairs_all(const std::string& simdir, const std::string& treedir,
    const std::string& envdir, const std::string& writepath,
    const snapnum_type snapnum_first, const snapnum_type snapnum_last) {

  // Parameters
  double rmin =  5.0;  // kpc/h
  double rmax = 20.0;  // kpc/h
  double velocity_threshold = 500.0;  // km/s
  float mstar_min = 0.704;  // 10^10 Msun/h

  // Get box size in ckpc/h; note conversion to float
  std::stringstream tmp_stream;
  tmp_stream << simdir << "/output/snapdir_" <<
      std::setfill('0') << std::setw(3) << snapnum_last << "/snap_" <<
      std::setfill('0') << std::setw(3) << snapnum_last << ".0.hdf5";
  std::string file_name = tmp_stream.str();
  real_type box_size = static_cast<real_type>(
      arepo::get_scalar_attribute<double>(file_name, "BoxSize"));

  // Get redshift of each snapshot.
  auto redshifts_all = cosmo::get_redshifts();

  // Read all overdensities and store them in a vector of vectors.
  std::cout << "Reading overdensities...\n";
  WallClock wall_clock;
  std::vector<std::vector<real_type>> overdensities(snapnum_last+1);
  for (auto snapnum = snapnum_first; snapnum <= snapnum_last; ++snapnum) {
    // Only proceed if there is at least one subhalo
    auto nsubs = subfind::get_scalar_attribute<uint32_t>(
        simdir + "/output", snapnum, "Nsubgroups_Total");
    if (nsubs == 0)
      continue;

    // Create filename
    tmp_stream.str("");
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
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";

//  // Load merger tree
//  wall_clock.start();
//  std::cout << "Loading merger tree...\n";
//  int filenum = -1;  // "concatenated" tree file
//  std::string name = "tree_extended";  // Full format
//  Tree tree(treedir, name, filenum);
//  std::cout << "Loaded merger tree. Total time: " <<
//      wall_clock.seconds() << " s.\n";
//  std::cout << "\n";

  // Iterate over snapshots
  std::cout << "Iterating over snapshots..." << std::endl;
  for (auto snapnum = snapnum_first; snapnum <= snapnum_last; ++snapnum) {
//    auto snap = tree.snapshot(snapnum);

    // To comoving coordinates
    real_type z = redshifts_all[snapnum];
    real_type rmin_comoving = rmin * (1.0 + z);
    real_type rmax_comoving = rmax * (1.0 + z);

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

//    // Iterate over subhalos and count mergers.
//    for (auto sub_it = snap.begin(); sub_it != snap.end(); ++sub_it) {
//      count_mergers_sub(*sub_it, writefile, overdensities);
//    }


    // Load galaxy positions and stellar masses.
    std::cout << "Loading data...\n";
    WallClock wall_clock;
    std::string basedir = simdir + "/output";
    auto sub_pos = subfind::read_block<Point>(
        basedir, snapnum, "Subhalo", "SubhaloPos", -1);
    auto sub_vel = subfind::read_block<Point>(
        basedir, snapnum, "Subhalo", "SubhaloVel", -1);
    auto sub_mstar = subfind::read_block<real_type>(
        basedir, snapnum, "Subhalo", "SubhaloMassType", parttype_stars);
    uint32_t nsubs = sub_pos.size();
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";

    // Create Galaxy objects and add to AllGalaxies container.
    std::cout << "Creating galaxy objects...\n";
    wall_clock.start();
    AllGalaxiesType ag;
    for (uint32_t sub_index = 0; sub_index < nsubs; ++sub_index) {
      // Add galaxy
      ag.add_galaxy(sub_pos[sub_index], GalaxyData(sub_index, sub_vel[sub_index]));

      // Add one "ghost" galaxy for each dimension (optimization:
      // only include ghost galaxies within 100 kpc/h of each border)
      // (Another possible optimization: implement 2D spatial search,
      // but probably none of this is necessary because of slow reading speed).
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

    // Iterate over real (non-ghost) galaxies
    std::cout << "Iterating over galaxies...\n";
    wall_clock.start();
    for (uint32_t i = 0; i < nsubs; ++i) {

      // Check that current galaxy is above minimum mass
      if (sub_mstar[i] < mstar_min)
        continue;

//      auto center_sub = tree.subhalo(snapnum, i);
      auto center_gal = ag.galaxy(i);

      auto neighbors = annulus_query(s, center_gal, rmin_comoving, rmax_comoving);

      for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
        auto other_gal = *it;
        if (within_velocity_threshold(center_gal, other_gal, velocity_threshold, z)) {
          // Do something, e.g., check mass ratio
          std::cout << "Inside the annulus!\n";
        }
      }

      // Print some output occasionally
      if ((i+1) % 10000 == 0) {
        std::cout << "Already processed " << i+1 << " galaxies.\n";
      }
    }
    std::cout << "Time: " << wall_clock.seconds() << " s.\n";





    // Flush and close file
    writefile.flush();
    writefile.close();
    std::cout << "Finished for snapshot " << snapnum << std::endl;
  }
}


int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 7) {
    std::cerr << "Usage: " << argv[0] << " simdir treedir envdir writepath" <<
        " snapnum_first snapnum_last\n";
    exit(1);
  }

  // Read input
  std::string simdir(argv[1]);
  std::string treedir(argv[2]);
  std::string envdir(argv[3]);
  std::string writepath(argv[4]);
  snapnum_type snapnum_first = atoi(argv[5]);
  snapnum_type snapnum_last = atoi(argv[6]);

  // Measure CPU and wall clock (real) time
  WallClock wall_clock;
  CPUClock cpu_clock;

  // Do stuff
  count_pairs_all(simdir, treedir, envdir, writepath, snapnum_first, snapnum_last);

  // Print wall clock time and speedup
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  std::cout << "Speedup: " << cpu_clock.seconds()/wall_clock.seconds() << ".\n";

  return 0;
}
