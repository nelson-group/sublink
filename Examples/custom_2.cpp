/** @file custom_2.cpp
 * @brief For a few galaxies, find all the times they underwent a merger
 *        with mass ratio > 1:4, and print the relevant Subfind IDs.
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
#include "../Util/SnapshotUtil.hpp"  // get_redshifts, get_times_Gyr
#include "../Util/TreeUtil.hpp"

static constexpr float major_merger_ratio = 1.0/4.0;

/** @brief Find the last major/minor mergers for a given subhalo. */
void merger_history_sub(Subhalo sub, std::ofstream& outfile) {

  // Print Subfind ID.
  outfile << "\n" << sub.index() << "\n";

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
      if ((merger_ratio >= major_merger_ratio) &&
          (merger_ratio <= 1.0/major_merger_ratio)) {
        // Print info.
        outfile << sub.snapnum() << "," << sub.index() << "," <<
            stmax_pair.first.snapnum() << "," << stmax_pair.first.index() << "," <<
            stmax_pair.second.snapnum() << "," << stmax_pair.second.index() << "\n";
      }
    }

    // Next iteration
    sub = first_prog;
    first_prog = sub.first_progenitor();
  }
}

void brendan_custom(const Tree& tree,
    const std::string& filebase) {

  // Read selected subhalos from first file
  std::ifstream infile(filebase + ".txt");
  std::vector<index_type> subfind_ids;
  index_type subf_id;
  while (infile >> subf_id) {
    subfind_ids.push_back(subf_id);
  }
  infile.close();
  snapnum_type snapnum = 135;

  // Iterate over subhalos
  std::ofstream outfile(filebase + "_output.txt");
  for (auto sub_index: subfind_ids) {
    auto sub = tree.subhalo(snapnum, sub_index);
    merger_history_sub(sub, outfile);
  }
  outfile.close();
}


int main()
{
  std::string treedir = "/n/ghernquist/vrodrigu/MergerTrees/output/Galaxies/Illustris/L75n1820FP";

  // Measure CPU and wall clock (real) time
  WallClock wall_clock_all;
  CPUClock cpu_clock_all;

  // Load merger tree
  WallClock wall_clock;
  std::cout << "Loading merger tree...\n";
  int filenum = -1;  // read from all merger tree files
  std::string name = "tree_extended";  // Full format
  Tree tree(treedir, name, filenum);
  std::cout << "Loaded merger tree. Total time: " <<
      wall_clock.seconds() << " s.\n";
  std::cout << "\n";

  // Do stuff for both files
  brendan_custom(tree, "isolated_ids");
  brendan_custom(tree, "pair_ids");

  // Print wall clock time and speedup
  std::cout << "Time: " << wall_clock_all.seconds() << " s.\n";
  std::cout << "Speedup: " << cpu_clock_all.seconds()/wall_clock.seconds() << ".\n";

  return 0;
}
