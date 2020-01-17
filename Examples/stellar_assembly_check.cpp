/** @file stellar_assembly_check.cpp
 * @brief For each galaxy from a given snapshot, estimate the ex situ fraction
 *        using the merger trees.
 *
 * @author Vicente Rodriguez-Gomez (v.rodriguez@irya.unam.mx)
 */

// Include some extra quantities from the merger trees:
#define COUNT_MERGERS

#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "../InputOutput/ReadSubfindHDF5.hpp"
#include "../InputOutput/ReadTreeHDF5.hpp"
#include "../InputOutput/GeneralHDF5.hpp"
#include "../Util/GeneralUtil.hpp"
#include "../Util/TreeUtil.hpp"

struct AccretedData {
  std::vector<real_type> StellarMassInstant;
  std::vector<real_type> StellarMassMax;
  std::vector<real_type> StellarMassInfall;

  /** Constructor.
   * @param[in] nsubs The number of subhalos in this snapshot.
   */
  AccretedData(const uint32_t nsubs)
      : StellarMassInstant(std::vector<real_type>(nsubs, 0)),
        StellarMassMax(std::vector<real_type>(nsubs, 0)),
        StellarMassInfall(std::vector<real_type>(nsubs, 0)) {
  }
};

/** @brief Get infall mass, considering all progenitors recursively.
 * Based on get_infall_pair function from TreeUtil.hpp.
 */
real_type get_infall_mass_recursive(Subhalo primary,
    Subhalo secondary) {
  assert(primary.is_valid() && secondary.is_valid());

  real_type running_total = 0;

  // Increment once just to start
  secondary = secondary.first_progenitor();

  while (true) {
    // If secondary branch is truncated, return.
    if (!secondary.is_valid())
      return running_total;

    // If primary branch is truncated, return.
    auto cur_snapnum = secondary.snapnum();
    primary = back_in_time(primary, cur_snapnum);
    if (!primary.is_valid())
      return running_total;

    // If primary skipped this snapshot, or is for some other reason found
    // at an earlier snapshot, "increment" secondary and try again.
    if (primary.snapnum() != cur_snapnum) {
      secondary = secondary.first_progenitor();
      continue;
    }

    // If we got here, then both primary and secondary are valid Subhalos at
    // the same snapshot.

    // Iterate over next progenitor
    for (auto next_prog = secondary.next_progenitor(); next_prog.is_valid();
        next_prog = next_prog.next_progenitor()) {
      running_total += get_infall_mass_recursive(primary, next_prog);
    }

    // If primary and secondary do not belong to the same FoF group,
    // we have reached infall.
    if (  primary.first_subhalo_in_fof_group() !=
        secondary.first_subhalo_in_fof_group()) {
      running_total += secondary.data().SubhaloMassType[4];
      return running_total;
    }

    // Increment secondary
    secondary = secondary.first_progenitor();
  }
  assert(false);
}

/** @brief Calculate ex situ fraction for a given subhalo. */
void check_sub(Subhalo sub, AccretedData& ad) {

  uint32_t index_orig = sub.index();

  // Iterate over first progenitor
  auto first_prog = sub.first_progenitor();
  while (first_prog.is_valid()) {

    // Iterate over next progenitor
    for (auto next_prog = first_prog.next_progenitor(); next_prog.is_valid();
        next_prog = next_prog.next_progenitor()) {

      // ------------------- INSTANTANEOUS -------------------------

      // Add "instantaneous" mass
      ad.StellarMassInstant[index_orig] += next_prog.data().SubhaloMassType[4];

      // ---------------------- STMAX ---------------------------

      // Progenitor properties at stellar tmax
      auto stmax_pair = get_stmax_pair(first_prog, next_prog);

      // Only proceed if stellar tmax is well defined
      // (NOTE: it *could* be possible that infall is well defined
      //  while tmax is not)
      auto snapnum_stmax = stmax_pair.second.snapnum();
      if (snapnum_stmax == -1)
        continue;

      // Only proceed if both stellar masses at stmax are > 0
      if (!((stmax_pair.first.data().SubhaloMassType[4] > 0) &&
          (stmax_pair.second.data().SubhaloMassType[4] > 0))) {
        continue;
      }

      // Add mass at stellar tmax
      real_type mass_secondary = std::min(
          stmax_pair.first.data().SubhaloMassType[4], stmax_pair.second.data().SubhaloMassType[4]);
      ad.StellarMassMax[index_orig] += mass_secondary;

      // ---------------------- INFALL --------------------------

//      // Progenitor properties at infall
//      auto infall_pair = get_infall_pair(first_prog, next_prog);
//
//      // Only proceed if infall is well defined
//      auto snapnum_infall = infall_pair.second.snapnum();
//      if (snapnum_infall == -1)
//        continue;
//
//      // Add mass at infall
//      ad.StellarMassInfall[index_orig] += infall_pair.second.data().SubhaloMassType[4];

      // Add masses at infall recursively.
      ad.StellarMassInfall[index_orig] += get_infall_mass_recursive(first_prog, next_prog);
    }

    // Next iteration
    sub = first_prog;
    first_prog = sub.first_progenitor();
  }
}

/** @brief Estimate ex situ stellar mass fractions. */
void stellar_assembly_check(const std::string& basedir, const std::string& treedir,
    const std::string& writepath, const snapnum_type snapnum) {

  // Load merger tree
  WallClock wall_clock;
  std::cout << "Loading merger tree...\n";
  int filenum = -1;  // read from all merger tree files
  std::string name = "tree_extended";  // Full format
  Tree tree(treedir, name, filenum);
  std::cout << "Loaded merger tree. Total time: " <<
      wall_clock.seconds() << " s.\n";
  std::cout << "\n";

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
    return;
  }

  // Store data here:
  AccretedData ad(nsubs);

  // Iterate over subhalos to calculate ex situ fraction.
  auto snap = tree.snapshot(snapnum);
  for (auto sub_it = snap.begin(); sub_it != snap.end(); ++sub_it) {
    check_sub(*sub_it, ad);
  }

  // Write to file.
  wall_clock.start();

  add_array(writefile, ad.StellarMassInstant, "StellarMassInstant", H5::PredType::NATIVE_FLOAT);
  add_array(writefile, ad.StellarMassMax, "StellarMassMax", H5::PredType::NATIVE_FLOAT);
  add_array(writefile, ad.StellarMassInfall, "StellarMassInfall", H5::PredType::NATIVE_FLOAT);

  // Close (and flush) file
  writefile.close();
  std::cout << "Finished for snapshot " << snapnum << std::endl;

}


int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0] << " basedir treedir writepath" <<
        " snapnum\n";
    exit(1);
  }

  // Read input
  std::string basedir(argv[1]);
  std::string treedir(argv[2]);
  std::string writepath(argv[3]);
  snapnum_type snapnum = atoi(argv[4]);

  // Measure CPU and wall clock (real) time
  WallClock wall_clock;
  CPUClock cpu_clock;

  // Do stuff
  stellar_assembly_check(basedir, treedir, writepath, snapnum);

  // Print wall clock time and speedup
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  std::cout << "Speedup: " << cpu_clock.seconds()/wall_clock.seconds() << ".\n";

  return 0;
}
