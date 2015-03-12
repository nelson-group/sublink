/** @file shy_custom_1.cpp
 * @brief For a few galaxies, find all the times they underwent a merger
 *        with mass ratio > 1:10, and the mass ratios of those mergers.
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
#include "../Util/SnapshotUtil.hpp"  // get_redshifts, get_times_Gyr
#include "../Util/TreeUtil.hpp"

static constexpr float minor_merger_ratio = 1.0/10.0;

/** @brief Find the last major/minor mergers for a given subhalo. */
void merger_history_sub(Subhalo sub) {

  // Print Subfind ID.
  std::cout << "\n" << sub.index() << "\n";

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

      // Check if minor merger
      float merger_ratio = stmax_pair.second.data().SubhaloMassType[4] / stmax_pair.first.data().SubhaloMassType[4];
      if ((merger_ratio >= minor_merger_ratio) &&
          (merger_ratio <= 1.0/minor_merger_ratio)) {
        // Print info.
        std::cout << sub.snapnum() << "," << std::min(merger_ratio, 1.0f/merger_ratio) << "\n";
      }
    }

    // Next iteration
    sub = first_prog;
    first_prog = sub.first_progenitor();
  }
}

/** @brief Get snapshot of last major (minor) merger. */
void shy_custom(const std::string& treedir) {

  // Selected subhalos:
  std::vector<index_type> subfind_ids = {219708, 262030, 377255, 267605, 195486, 386640};
  snapnum_type snapnum = 135;

  // Load merger tree
  WallClock wall_clock;
  std::cout << "Loading merger tree...\n";
  int filenum = -1;  // "concatenated" tree file
  std::string name = "tree_extended";  // Full format
  Tree tree(treedir, name, filenum);
  std::cout << "Loaded merger tree. Total time: " <<
      wall_clock.seconds() << " s.\n";
  std::cout << "\n";

  // Iterate over subhalos
  for (auto sub_index: subfind_ids) {
    auto sub = tree.subhalo(snapnum, sub_index);
    merger_history_sub(sub);
  }
}


int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " treedir\n";
    exit(1);
  }

  // Read input
  std::string treedir(argv[1]);

  // Measure CPU and wall clock (real) time
  WallClock wall_clock;
  CPUClock cpu_clock;

  // Do stuff
  shy_custom(treedir);

  // Print wall clock time and speedup
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  std::cout << "Speedup: " << cpu_clock.seconds()/wall_clock.seconds() << ".\n";

  return 0;
}
