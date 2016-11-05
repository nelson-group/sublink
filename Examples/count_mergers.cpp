/** @file count_mergers.cpp
 * @brief Count mergers and print their main properties into files.
 *
 * @author Vicente Rodriguez-Gomez (vrg@jhu.edu)
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

/** @brief Count mergers for a given subhalo and print data to file.
 * @param[in] sub The subhalo of interest.
 * @param[in] writefile An open file to write the results.
 */
void count_mergers_sub(Subhalo sub, std::ofstream& writefile,
    const std::vector<std::vector<real_type>>& overdensities) {
  // Only proceed if @a sub has at least one progenitor.
  auto first_prog = sub.first_progenitor();
  if (!first_prog.is_valid())
    return;

  // Iterate over "next progenitor" link to count mergers.
  for (auto next_prog = first_prog.next_progenitor();
      next_prog.is_valid(); next_prog = next_prog.next_progenitor()) {

    // Descendant properties.
    auto is_central = static_cast<int>(sub.data().FirstSubhaloInFOFGroupID == sub.data().SubhaloID);
    writefile << std::setprecision(10) <<
        sub.data().SubhaloMassType[4] << "," <<
        sub.data().Mass << "," <<
        sub.data().SubhaloSFR << "," <<
        sub.data().SubhaloStellarPhotometrics[4] - sub.data().SubhaloStellarPhotometrics[5] << "," <<
        overdensities[sub.data().SnapNum][sub.data().SubfindID] << "," <<
        is_central << ",";

    // Progenitor properties right before merger.
    writefile << std::setprecision(10) <<
        first_prog.data().SubhaloMassType[4] << "," <<
        first_prog.data().Mass << "," <<
        first_prog.data().SubhaloSFR << "," <<
        first_prog.data().SubhaloStellarPhotometrics[4] - first_prog.data().SubhaloStellarPhotometrics[5] << "," <<
        overdensities[first_prog.data().SnapNum][first_prog.data().SubfindID] << "," <<
        next_prog.data().SubhaloMassType[4] << "," <<
        next_prog.data().Mass << "," <<
        next_prog.data().SubhaloSFR << "," <<
        next_prog.data().SubhaloStellarPhotometrics[4] - next_prog.data().SubhaloStellarPhotometrics[5] << "," <<
        overdensities[next_prog.data().SnapNum][next_prog.data().SubfindID] << ",";

    // Progenitor properties at stellar tmax.
    real_type mstar_stmax_1 = -1;
    real_type mgal_stmax_1 = -1;
    real_type sfr_stmax_1 = -1;
    real_type g_minus_r_stmax_1 = -1;
    real_type overdensity_stmax_1 = -1;
    real_type mstar_stmax_2 = -1;
    real_type mgal_stmax_2 = -1;
    real_type sfr_stmax_2 = -1;
    real_type g_minus_r_stmax_2 = -1;
    real_type overdensity_stmax_2 = -1;
    snapnum_type snapnum_stmax = -1;
    auto stmax_pair = get_stmax_pair(first_prog, next_prog);
    if (stmax_pair.first.is_valid()) {
      mstar_stmax_1 = stmax_pair.first.data().SubhaloMassType[4];
      mgal_stmax_1 = stmax_pair.first.data().Mass;
      sfr_stmax_1 = stmax_pair.first.data().SubhaloSFR;
      g_minus_r_stmax_1 = stmax_pair.first.data().SubhaloStellarPhotometrics[4] - stmax_pair.first.data().SubhaloStellarPhotometrics[5];
      overdensity_stmax_1 = overdensities[stmax_pair.first.data().SnapNum][stmax_pair.first.data().SubfindID];
      mstar_stmax_2 = stmax_pair.second.data().SubhaloMassType[4];
      mgal_stmax_2 = stmax_pair.second.data().Mass;
      sfr_stmax_2 = stmax_pair.second.data().SubhaloSFR;
      g_minus_r_stmax_2 = stmax_pair.second.data().SubhaloStellarPhotometrics[4] - stmax_pair.second.data().SubhaloStellarPhotometrics[5];
      overdensity_stmax_2 = overdensities[stmax_pair.second.data().SnapNum][stmax_pair.second.data().SubfindID];
      snapnum_stmax = stmax_pair.second.data().SnapNum;
    }
    writefile << std::setprecision(10) <<
        mstar_stmax_1 << "," <<
        mgal_stmax_1 << "," <<
        sfr_stmax_1 << "," <<
        g_minus_r_stmax_1 << "," <<
        overdensity_stmax_1 << "," <<
        mstar_stmax_2 << "," <<
        mgal_stmax_2 << "," <<
        sfr_stmax_2 << "," <<
        g_minus_r_stmax_2 << "," <<
        overdensity_stmax_2 << "," <<
        snapnum_stmax << ",";

    // Progenitor properties at infall.
    real_type mstar_infall_1 = -1;
    real_type mgal_infall_1 = -1;
    real_type sfr_infall_1 = -1;
    real_type g_minus_r_infall_1 = -1;
    real_type overdensity_infall_1 = -1;
    real_type mstar_infall_2 = -1;
    real_type mgal_infall_2 = -1;
    real_type sfr_infall_2 = -1;
    real_type g_minus_r_infall_2 = -1;
    real_type overdensity_infall_2 = -1;
    snapnum_type snapnum_infall = -1;
    auto infall_pair = get_infall_pair(first_prog, next_prog);
    if (infall_pair.first.is_valid()) {
      mstar_infall_1 = infall_pair.first.data().SubhaloMassType[4];
      mgal_infall_1 = infall_pair.first.data().Mass;
      sfr_infall_1 = infall_pair.first.data().SubhaloSFR;
      g_minus_r_infall_1 = infall_pair.first.data().SubhaloStellarPhotometrics[4] - infall_pair.first.data().SubhaloStellarPhotometrics[5];
      overdensity_infall_1 = overdensities[infall_pair.first.data().SnapNum][infall_pair.first.data().SubfindID];
      mstar_infall_2 = infall_pair.second.data().SubhaloMassType[4];
      mgal_infall_2 = infall_pair.second.data().Mass;
      sfr_infall_2 = infall_pair.second.data().SubhaloSFR;
      g_minus_r_infall_2 = infall_pair.second.data().SubhaloStellarPhotometrics[4] - infall_pair.second.data().SubhaloStellarPhotometrics[5];
      overdensity_infall_2 = overdensities[infall_pair.second.data().SnapNum][infall_pair.second.data().SubfindID];
      snapnum_infall = infall_pair.second.data().SnapNum;
    }
    writefile << std::setprecision(10) <<
        mstar_infall_1 << "," <<
        mgal_infall_1 << "," <<
        sfr_infall_1 << "," <<
        g_minus_r_infall_1 << "," <<
        overdensity_infall_1 << "," <<
        mstar_infall_2 << "," <<
        mgal_infall_2 << "," <<
        sfr_infall_2 << "," <<
        g_minus_r_infall_2 << "," <<
        overdensity_infall_2 << "," <<
        snapnum_infall << "\n";
  }
}

/** @brief Count mergers and print to files. */
void count_mergers_all(const std::string& simdir, const std::string& treedir,
    const std::string& writepath, const snapnum_type snapnum_first,
    const snapnum_type snapnum_last) {

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
    std::stringstream tmp_stream;
    tmp_stream << simdir << "/postprocessing/environment/environment_" <<
        std::setfill('0') << std::setw(3) << snapnum << ".hdf5";
    std::string filename = tmp_stream.str();

    // Quickly check if file exists
    std::ifstream file(filename);
    if (!file) {
      // No overdensities to read. Fill array with -1001's
      overdensities[snapnum].resize(nsubs, -1001);
      continue;
    }
    file.close();

    // Read overdensities
    overdensities[snapnum] = read_dataset<real_type>(filename, "delta");
    assert(overdensities[snapnum].size() == nsubs);
  }
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";

  // Load merger tree
  wall_clock.start();
  std::cout << "Loading merger tree...\n";
  int filenum = -1;  // read from all merger tree files
  std::string name = "tree_extended";  // Full format
  Tree tree(treedir, name, filenum);
  std::cout << "Loaded merger tree. Total time: " <<
      wall_clock.seconds() << " s.\n";
  std::cout << "\n";

  // Iterate over snapshots
  std::cout << "Iterating over snapshots..." << std::endl;
  for (auto snapnum = snapnum_first; snapnum <= snapnum_last; ++snapnum) {
    auto snap = tree.snapshot(snapnum);

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

    // Iterate over subhalos and count mergers.
    for (auto sub_it = snap.begin(); sub_it != snap.end(); ++sub_it) {
      count_mergers_sub(*sub_it, writefile, overdensities);
    }

    // Flush and close file
    writefile.flush();
    writefile.close();
    std::cout << "Finished for snapshot " << snapnum << std::endl;
  }
}


int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 6) {
    std::cerr << "Usage: " << argv[0] << " simdir treedir writepath" <<
        " snapnum_first snapnum_last\n";
    exit(1);
  }

  // Read input
  std::string simdir(argv[1]);
  std::string treedir(argv[2]);
  std::string writepath(argv[3]);
  snapnum_type snapnum_first = atoi(argv[4]);
  snapnum_type snapnum_last = atoi(argv[5]);

  // Measure CPU and wall clock (real) time
  WallClock wall_clock;
  CPUClock cpu_clock;

  // Do stuff
  count_mergers_all(simdir, treedir, writepath, snapnum_first, snapnum_last);

  // Print wall clock time and speedup
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  std::cout << "Speedup: " << cpu_clock.seconds()/wall_clock.seconds() << ".\n";

  return 0;
}
