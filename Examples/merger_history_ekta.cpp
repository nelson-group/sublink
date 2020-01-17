/** @file merger_history.cpp
 * @brief For each galaxy, print some info about the last merger.
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
#include "../Util/Cosmology.hpp"  // cosmo::get_redshifts, cosmo::get_times_Gyr
#include "../Util/TreeUtil.hpp"

static constexpr float major_merger_ratio = 1.0/4.0;
static constexpr float minor_merger_ratio = 1.0/10.0;

/** @brief Type for number of mergers. */
typedef uint32_t count_type;

struct MergerData {

  // Major mergers
  std::vector<snapnum_type> SnapNumLastMajorMerger;
  std::vector<snapnum_type> SnapNumBeforeLastMajorMerger;
  std::vector<index_type> SubfindIDLastMajorMerger;

  // Minor mergers
  std::vector<snapnum_type> SnapNumLastMinorMerger;
  std::vector<snapnum_type> SnapNumBeforeLastMinorMerger;
  std::vector<index_type> SubfindIDLastMinorMerger;

  // All mergers
  std::vector<snapnum_type> SnapNumLastMerger;
  std::vector<snapnum_type> SnapNumBeforeLastMerger;
  std::vector<index_type> SubfindIDLastMerger;

  /** Constructor.
   * @param[in] nsubs The number of subhalos in this snapshot.
   */
  MergerData(const uint32_t nsubs)
      : SnapNumLastMajorMerger(std::vector<snapnum_type>(nsubs, -1)),
        SnapNumBeforeLastMajorMerger(std::vector<snapnum_type>(nsubs, -1)),
        SubfindIDLastMajorMerger(std::vector<index_type>(nsubs, -1)),

        SnapNumLastMinorMerger(std::vector<snapnum_type>(nsubs, -1)),
        SnapNumBeforeLastMinorMerger(std::vector<snapnum_type>(nsubs, -1)),
        SubfindIDLastMinorMerger(std::vector<index_type>(nsubs, -1)),

        SnapNumLastMerger(std::vector<snapnum_type>(nsubs, -1)),
        SnapNumBeforeLastMerger(std::vector<snapnum_type>(nsubs, -1)),
        SubfindIDLastMerger(std::vector<index_type>(nsubs, -1)) {
  }
};

/** @brief Find the last major/minor mergers for a given subhalo. */
void merger_history_sub(Subhalo sub, MergerData& md) {

  uint32_t index_orig = sub.index();

  bool no_major_merger_yet = true;
  bool no_minor_merger_yet = true;
  bool no_merger_yet = true;

  // Iterate over first progenitor
  auto first_prog = sub.first_progenitor();
  while (first_prog.is_valid()) {

    // Iterate over next progenitor
    for (auto next_prog = first_prog.next_progenitor(); next_prog.is_valid();
        next_prog = next_prog.next_progenitor()) {

      // Progenitor properties at stellar tmax
      auto stmax_pair = get_stmax_pair(first_prog, next_prog);

      // Progenitor properties at infall
      auto infall_pair = get_infall_pair(first_prog, next_prog);

      // Only proceed if infall is well defined
      auto snapnum_infall = infall_pair.second.snapnum();
      if (snapnum_infall == -1)
        continue;

      // Only proceed if stellar tmax is well defined
      auto snapnum_stmax = stmax_pair.second.snapnum();
      if (snapnum_stmax == -1)
        continue;

      // Only proceed if both stellar masses at stmax are > 0
      if (!((stmax_pair.first.data().SubhaloMassType[4] > 0) &&
          (stmax_pair.second.data().SubhaloMassType[4] > 0))) {
        continue;
      }

      // Calculate stellar mass ratio (can be > 1)
      real_type mass_ratio = stmax_pair.second.data().SubhaloMassType[4] / stmax_pair.first.data().SubhaloMassType[4];

      // Check if major merger
      if ((mass_ratio >= major_merger_ratio) &&
          (mass_ratio <= 1.0/major_merger_ratio)) {
        // Check if most recent
        if (no_major_merger_yet) {
          md.SnapNumLastMajorMerger[index_orig] = sub.snapnum();
          md.SnapNumBeforeLastMajorMerger[index_orig] = next_prog.snapnum();
          md.SubfindIDLastMajorMerger[index_orig] = next_prog.index();
          no_major_merger_yet = false;
        }
      }

      // Check if minor merger (does not overlap with major mergers)
      if (((mass_ratio >= minor_merger_ratio) &&
           (mass_ratio <  major_merger_ratio)) ||
          ((mass_ratio >  1.0/major_merger_ratio) &&
           (mass_ratio <= 1.0/minor_merger_ratio))) {
        // Check if most recent
        if (no_minor_merger_yet) {
          md.SnapNumLastMinorMerger[index_orig] = sub.snapnum();
          md.SnapNumBeforeLastMinorMerger[index_orig] = next_prog.snapnum();
          md.SubfindIDLastMinorMerger[index_orig] = next_prog.index();
          no_minor_merger_yet = false;
        }
      }

      // All mergers
      if (true) {
        // Check if most recent
        if (no_merger_yet) {
          md.SnapNumLastMerger[index_orig] = sub.snapnum();
          md.SnapNumBeforeLastMerger[index_orig] = next_prog.snapnum();
          md.SubfindIDLastMerger[index_orig] = next_prog.index();
          no_merger_yet = false;
        }
      }
    }

    // Next iteration
    sub = first_prog;
    first_prog = sub.first_progenitor();
  }
}

/** @brief Get snapshot of last major (minor) merger. */
void merger_history_all(
    const std::string& suite, const std::string& basedir,
    const std::string& treedir, const std::string& writepath,
    const snapnum_type snapnum_first, const snapnum_type snapnum_last) {

  (void) suite;  // Silence unused variable warning

  // Load merger tree
  WallClock wall_clock;
  std::cout << "Loading merger tree...\n";
  int filenum = -1;  // read from all merger tree files
  std::string name = "tree_extended";  // Full format
  Tree tree(treedir, name, filenum);
  std::cout << "Loaded merger tree. Total time: " <<
      wall_clock.seconds() << " s.\n";
  std::cout << "\n";

  // Iterate over snapshots (backwards)
  wall_clock.start();
  std::cout << "Iterating over snapshots..." << std::endl;
  for (auto snapnum = snapnum_last; snapnum >= snapnum_first; --snapnum) {

    // AD HOC: if L75n1820FP, skip snapshots 53 and 55.
    if (basedir == "/n/hernquistfs1/Illustris/Runs/L75n1820FP/output") {
      if ((snapnum == 53) || (snapnum == 55)) {
        std::cout << "WARNING: Skipping snapshot " << snapnum << ".\n";
        continue;
      }
    }

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

    // Store data here:
    MergerData md(nsubs);

    // Iterate over subhalos and find the last major (minor) merger.
    auto snap = tree.snapshot(snapnum);
    for (auto sub_it = snap.begin(); sub_it != snap.end(); ++sub_it) {
      merger_history_sub(*sub_it, md);
    }

    // Write to file.
    wall_clock.start();

    add_array(writefile, md.SnapNumLastMajorMerger, "SnapNumLastMajorMerger", H5::PredType::NATIVE_INT16);
    add_array(writefile, md.SnapNumBeforeLastMajorMerger, "SnapNumBeforeLastMajorMerger", H5::PredType::NATIVE_INT16);
    add_array(writefile, md.SubfindIDLastMajorMerger, "SubfindIDLastMajorMerger", H5::PredType::NATIVE_INT32);

    add_array(writefile, md.SnapNumLastMinorMerger, "SnapNumLastMinorMerger", H5::PredType::NATIVE_INT16);
    add_array(writefile, md.SnapNumBeforeLastMinorMerger, "SnapNumBeforeLastMinorMerger", H5::PredType::NATIVE_INT16);
    add_array(writefile, md.SubfindIDLastMinorMerger, "SubfindIDLastMinorMerger", H5::PredType::NATIVE_INT32);

    add_array(writefile, md.SnapNumLastMerger, "SnapNumLastMerger", H5::PredType::NATIVE_INT16);
    add_array(writefile, md.SnapNumBeforeLastMerger, "SnapNumBeforeLastMerger", H5::PredType::NATIVE_INT16);
    add_array(writefile, md.SubfindIDLastMerger, "SubfindIDLastMerger", H5::PredType::NATIVE_INT32);

    // Close (and flush) file
    writefile.close();
    std::cout << "Finished for snapshot " << snapnum << std::endl;
  }
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";
}


int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 7) {
    std::cerr << "Usage: " << argv[0] << " suite basedir treedir writepath" <<
        " snapnum_first snapnum_last\n";
    exit(1);
  }

  // Read input
  std::string suite(argv[1]);
  std::string basedir(argv[2]);
  std::string treedir(argv[3]);
  std::string writepath(argv[4]);
  snapnum_type snapnum_first = atoi(argv[5]);
  snapnum_type snapnum_last = atoi(argv[6]);

  // Measure CPU and wall clock (real) time
  WallClock wall_clock;
  CPUClock cpu_clock;

  // Do stuff
  merger_history_all(suite, basedir, treedir, writepath, snapnum_first, snapnum_last);

  // Print wall clock time and CPU time
  std::cout << "Wall clock time: " << wall_clock.seconds() << " s.\n";
  std::cout << "CPU time: " << cpu_clock.seconds() << " s.\n";

  return 0;
}
