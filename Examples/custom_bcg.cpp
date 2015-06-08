/** @file custom_bcg.cpp
 * @brief Print some info for the most massive galaxies in the simulation.
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
static constexpr float major_merger_ratio = 1.0/4.0;

/** @brief Print the contributions from major/minor/all mergers for each snapshot. */
void merger_history_sub(Subhalo sub) {

  // Print Subfind ID.
  std::cout << sub.index() << "\n";

  // Iterate over first progenitor
  auto first_prog = sub.first_progenitor();
  while (first_prog.is_valid()) {

    // Print the following info for each step:
    // snapnum, mstar, sfr, mstar_from_major_mergers, mstar_from_minor_mergers,
    // mstar_from_all_mergers

    std::cout << sub.snapnum() << ",";
    std::cout << sub.data().SubhaloMassType[4] << ",";
    std::cout << sub.data().SubhaloSFR << ",";

    // Mass accretion from current snapshot
    real_type mstar_from_major_mergers = 0.0;
    real_type mstar_from_minor_mergers = 0.0;
    real_type mstar_from_mergers = 0.0;

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

      real_type mass_primary = std::max(
          stmax_pair.first.data().SubhaloMassType[4], stmax_pair.second.data().SubhaloMassType[4]);
      real_type mass_secondary = std::min(
          stmax_pair.first.data().SubhaloMassType[4], stmax_pair.second.data().SubhaloMassType[4]);

      // Check if major merger
      if (mass_secondary/mass_primary >= major_merger_ratio) {
        mstar_from_major_mergers += mass_secondary;
      }
      // Check if minor merger
      if (mass_secondary/mass_primary >= minor_merger_ratio) {
        mstar_from_minor_mergers += mass_secondary;
      }
      // Any merger
      mstar_from_mergers += mass_secondary;
    }

    std::cout << mstar_from_major_mergers << ",";
    std::cout << mstar_from_minor_mergers << ",";
    std::cout << mstar_from_mergers << "\n";

    // Next iteration
    sub = first_prog;
    first_prog = sub.first_progenitor();
  }

  // Finished for this subhalo.
  std::cout << "\n";

}

/** @brief Get snapshot of last major (minor) merger. */
void custom_bcg(const std::string& treedir) {

  snapnum_type snapnum = 135;

  // Selected subhalos (10 largest mstar at z=0)
  std::vector<index_type> subfind_ids = {
      0,  30430,  66080,  59384,  16937,  80734,  93165, 142714, 99148, 41088};

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

int main()
{
  // Illustris-1
  std::string treedir = "/n/ghernquist/vrodrigu/MergerTrees/output/Galaxies/Illustris/L75n1820FP";

  // Measure CPU and wall clock (real) time
  WallClock wall_clock;
  CPUClock cpu_clock;

  // Do stuff
  custom_bcg(treedir);

  // Print wall clock time and speedup
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  std::cout << "Speedup: " << cpu_clock.seconds()/wall_clock.seconds() << ".\n";

  return 0;
}
