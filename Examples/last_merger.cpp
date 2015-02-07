/** @file last_merger.cpp
 * @brief Determine the snapshot number when each galaxy had its
 *        last major or minor merger, for two different mass ratio
 *        definitions.
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

#include "../InputOutput/ReadSubfindHDF5.hpp"
#include "../InputOutput/ReadTreeHDF5.hpp"
#include "../InputOutput/GeneralHDF5.hpp"
#include "../Util/GeneralUtil.hpp"
#include "../Util/TreeUtil.hpp"

static constexpr float major_merger_ratio = 1.0/4.0;
static constexpr float minor_merger_ratio = 1.0/10.0;

/** @brief Find the last major/minor mergers for a given subhalo. */
void last_major_merger_sub(Subhalo sub,
    std::vector<snapnum_type>& SnapNumLastMajorMerger,
    std::vector<snapnum_type>& SnapNumLastMajorMergerBaryonic,
    std::vector<snapnum_type>& SnapNumLastMinorMerger,
    std::vector<snapnum_type>& SnapNumLastMinorMergerBaryonic) {

  uint32_t orig_index = sub.index();

  bool no_major_merger_yet = true;
  bool no_major_merger_yet_baryonic = true;
  bool no_minor_merger_yet = true;
  bool no_minor_merger_yet_baryonic = true;

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

      ///////////////////
      // GALAXY MASSES //
      ///////////////////

      // Only proceed if both baryonic masses at stmax are > 0
      // (note that baryonic mass = 0 implies stellar mass = 0)
      if (!((stmax_pair.first.data().Mass > 0) &&
          (stmax_pair.second.data().Mass > 0))) {
        continue;
      }

      if ((no_major_merger_yet_baryonic) &&
          (stmax_pair.second.data().Mass / stmax_pair.first.data().Mass >= major_merger_ratio) &&
          (stmax_pair.second.data().Mass / stmax_pair.first.data().Mass <= 1.0/major_merger_ratio)) {
        SnapNumLastMajorMergerBaryonic[orig_index] = sub.snapnum();
        no_major_merger_yet_baryonic = false;
      }

      if ((no_minor_merger_yet_baryonic) &&
          (stmax_pair.second.data().Mass / stmax_pair.first.data().Mass >= minor_merger_ratio) &&
          (stmax_pair.second.data().Mass / stmax_pair.first.data().Mass <= 1.0/minor_merger_ratio)) {
        SnapNumLastMinorMergerBaryonic[orig_index] = sub.snapnum();
        no_minor_merger_yet_baryonic = false;
      }

      ///////////////////
      // STELLAR MASSES //
      ///////////////////

      // Only proceed if both stellar masses at stmax are > 0
      if (!((stmax_pair.first.data().SubhaloMassType[4] > 0) &&
          (stmax_pair.second.data().SubhaloMassType[4] > 0))) {
        continue;
      }

      if ((no_major_merger_yet) &&
          (stmax_pair.second.data().SubhaloMassType[4] / stmax_pair.first.data().SubhaloMassType[4] >= major_merger_ratio) &&
          (stmax_pair.second.data().SubhaloMassType[4] / stmax_pair.first.data().SubhaloMassType[4] <= 1.0/major_merger_ratio)) {
        SnapNumLastMajorMerger[orig_index] = sub.snapnum();
        no_major_merger_yet = false;
      }

      if ((no_minor_merger_yet) &&
          (stmax_pair.second.data().SubhaloMassType[4] / stmax_pair.first.data().SubhaloMassType[4] >= minor_merger_ratio) &&
          (stmax_pair.second.data().SubhaloMassType[4] / stmax_pair.first.data().SubhaloMassType[4] <= 1.0/minor_merger_ratio)) {
        SnapNumLastMinorMerger[orig_index] = sub.snapnum();
        no_minor_merger_yet = false;
      }
    }

    // Next iteration
    sub = first_prog;
    first_prog = sub.first_progenitor();
  }
}

/** @brief Get snapshot of last major (minor) merger. */
void last_major_merger_all(const std::string& basedir, const std::string& treedir,
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

  // Iterate over snapshots
  wall_clock.start();
  std::cout << "Iterating over snapshots..." << std::endl;
  for (auto snapnum = snapnum_first; snapnum <= snapnum_last; ++snapnum) {

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

    // Store answer in these arrays, initialized to-1
    std::vector<snapnum_type> SnapNumLastMajorMerger(nsubs, -1);
    std::vector<snapnum_type> SnapNumLastMajorMergerBaryonic(nsubs, -1);
    std::vector<snapnum_type> SnapNumLastMinorMerger(nsubs, -1);
    std::vector<snapnum_type> SnapNumLastMinorMergerBaryonic(nsubs, -1);

    // Iterate over subhalos and find the last major (minor) merger.
    auto snap = tree.snapshot(snapnum);
    for (auto sub_it = snap.begin(); sub_it != snap.end(); ++sub_it) {
      last_major_merger_sub(*sub_it,
          SnapNumLastMajorMerger, SnapNumLastMajorMergerBaryonic,
          SnapNumLastMinorMerger, SnapNumLastMinorMergerBaryonic);
    }

    // Write to file.
    wall_clock.start();
    add_array(writefile, SnapNumLastMajorMerger, "SnapNumLastMajorMerger",
        H5::PredType::NATIVE_INT16);
    add_array(writefile, SnapNumLastMajorMergerBaryonic, "SnapNumLastMajorMergerBaryonic",
        H5::PredType::NATIVE_INT16);
    add_array(writefile, SnapNumLastMinorMerger, "SnapNumLastMinorMerger",
        H5::PredType::NATIVE_INT16);
    add_array(writefile, SnapNumLastMinorMergerBaryonic, "SnapNumLastMinorMergerBaryonic",
        H5::PredType::NATIVE_INT16);

    // Close (and flush) file
    writefile.close();
    std::cout << "Finished for snapshot " << snapnum << std::endl;
  }
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";
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
  snapnum_type snapnum_first = atoi(argv[4]);
  snapnum_type snapnum_last = atoi(argv[5]);

  // Measure CPU and wall clock (real) time
  WallClock wall_clock;
  CPUClock cpu_clock;

  // Do stuff
  last_major_merger_all(basedir, treedir, writepath, snapnum_first, snapnum_last);

  // Print wall clock time and speedup
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  std::cout << "Speedup: " << cpu_clock.seconds()/wall_clock.seconds() << ".\n";

  return 0;
}
