/** @file merger_history.cpp
 * @brief For each galaxy, print some simple merger statistics, such as
 *        the snapshot number of the last major merger and the number of
 *        major mergers in the last Gyr.
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

static constexpr float major_merger_ratio = 1.0/4.0;
static constexpr float minor_merger_ratio = 1.0/10.0;

/** @brief Type for number of mergers. */
typedef uint32_t count_type;

struct MergerData {
  //////////////////
  // STELLAR MASS //
  //////////////////

  // Major mergers
  std::vector<snapnum_type> SnapNumLastMajorMerger;
  std::vector<count_type> NumMajorMergersLastGyr;
  std::vector<count_type> NumMajorMergersTotal;
  std::vector<count_type> NumMajorMergersSinceRedshiftOne;
  std::vector<count_type> NumMajorMergersSinceRedshiftTwo;
  std::vector<count_type> NumMajorMergersSinceRedshiftThree;

  // Minor mergers
  std::vector<snapnum_type> SnapNumLastMinorMerger;
  std::vector<count_type> NumMinorMergersLastGyr;
  std::vector<count_type> NumMinorMergersTotal;
  std::vector<count_type> NumMinorMergersSinceRedshiftOne;
  std::vector<count_type> NumMinorMergersSinceRedshiftTwo;
  std::vector<count_type> NumMinorMergersSinceRedshiftThree;

  ///////////////////
  // BARYONIC MASS //
  ///////////////////

  // Major mergers
  std::vector<snapnum_type> SnapNumLastMajorMergerBaryonic;
  std::vector<count_type> NumMajorMergersLastGyrBaryonic;
  std::vector<count_type> NumMajorMergersTotalBaryonic;
  std::vector<count_type> NumMajorMergersSinceRedshiftOneBaryonic;
  std::vector<count_type> NumMajorMergersSinceRedshiftTwoBaryonic;
  std::vector<count_type> NumMajorMergersSinceRedshiftThreeBaryonic;

  // Minor mergers
  std::vector<snapnum_type> SnapNumLastMinorMergerBaryonic;
  std::vector<count_type> NumMinorMergersLastGyrBaryonic;
  std::vector<count_type> NumMinorMergersTotalBaryonic;
  std::vector<count_type> NumMinorMergersSinceRedshiftOneBaryonic;
  std::vector<count_type> NumMinorMergersSinceRedshiftTwoBaryonic;
  std::vector<count_type> NumMinorMergersSinceRedshiftThreeBaryonic;

  /** Constructor.
   * @param[in] nsubs The number of subhalos in this snapshot.
   */
  MergerData(const uint32_t nsubs)
      : SnapNumLastMajorMerger(std::vector<snapnum_type>(nsubs, -1)),
        NumMajorMergersLastGyr(std::vector<count_type>(nsubs, 0)),
        NumMajorMergersTotal(std::vector<count_type>(nsubs, 0)),
        NumMajorMergersSinceRedshiftOne(std::vector<count_type>(nsubs, 0)),
        NumMajorMergersSinceRedshiftTwo(std::vector<count_type>(nsubs, 0)),
        NumMajorMergersSinceRedshiftThree(std::vector<count_type>(nsubs, 0)),

        SnapNumLastMinorMerger(std::vector<snapnum_type>(nsubs, -1)),
        NumMinorMergersLastGyr(std::vector<count_type>(nsubs, 0)),
        NumMinorMergersTotal(std::vector<count_type>(nsubs, 0)),
        NumMinorMergersSinceRedshiftOne(std::vector<count_type>(nsubs, 0)),
        NumMinorMergersSinceRedshiftTwo(std::vector<count_type>(nsubs, 0)),
        NumMinorMergersSinceRedshiftThree(std::vector<count_type>(nsubs, 0)),

        SnapNumLastMajorMergerBaryonic(std::vector<snapnum_type>(nsubs, -1)),
        NumMajorMergersLastGyrBaryonic(std::vector<count_type>(nsubs, 0)),
        NumMajorMergersTotalBaryonic(std::vector<count_type>(nsubs, 0)),
        NumMajorMergersSinceRedshiftOneBaryonic(std::vector<count_type>(nsubs, 0)),
        NumMajorMergersSinceRedshiftTwoBaryonic(std::vector<count_type>(nsubs, 0)),
        NumMajorMergersSinceRedshiftThreeBaryonic(std::vector<count_type>(nsubs, 0)),

        SnapNumLastMinorMergerBaryonic(std::vector<snapnum_type>(nsubs, -1)),
        NumMinorMergersLastGyrBaryonic(std::vector<count_type>(nsubs, 0)),
        NumMinorMergersTotalBaryonic(std::vector<count_type>(nsubs, 0)),
        NumMinorMergersSinceRedshiftOneBaryonic(std::vector<count_type>(nsubs, 0)),
        NumMinorMergersSinceRedshiftTwoBaryonic(std::vector<count_type>(nsubs, 0)),
        NumMinorMergersSinceRedshiftThreeBaryonic(std::vector<count_type>(nsubs, 0)) {
  }
};

/** @brief Find the last major/minor mergers for a given subhalo. */
void merger_history_sub(Subhalo sub,
    const std::vector<real_type>& redshifts_all,
    const std::vector<real_type>& times_all,
    MergerData& md) {

  uint32_t index_orig = sub.index();
  snapnum_type snapnum_orig = sub.snapnum();

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

      /////////////////////
      // BARYONIC MASSES //
      /////////////////////

      // Only proceed if both baryonic (galaxy) masses at stmax are > 0
      // (note that baryonic mass = 0 implies stellar mass = 0)
      if (!((stmax_pair.first.data().Mass > 0) &&
          (stmax_pair.second.data().Mass > 0))) {
        continue;
      }

      // Check if major merger
      if ((stmax_pair.second.data().Mass / stmax_pair.first.data().Mass >= major_merger_ratio) &&
          (stmax_pair.second.data().Mass / stmax_pair.first.data().Mass <= 1.0/major_merger_ratio)) {
        // Check if most recent
        if (no_major_merger_yet_baryonic) {
          md.SnapNumLastMajorMergerBaryonic[index_orig] = sub.snapnum();
          no_major_merger_yet_baryonic = false;
        }
        // Add to merger counters
        md.NumMajorMergersTotalBaryonic[index_orig] += 1;
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 1.0)
          md.NumMajorMergersLastGyrBaryonic[index_orig] += 1;
        if (redshifts_all[sub.snapnum()] < 1.0)
          md.NumMajorMergersSinceRedshiftOneBaryonic[index_orig] += 1;
        if (redshifts_all[sub.snapnum()] < 2.0)
          md.NumMajorMergersSinceRedshiftTwoBaryonic[index_orig] += 1;
        if (redshifts_all[sub.snapnum()] < 3.0)
          md.NumMajorMergersSinceRedshiftThreeBaryonic[index_orig] += 1;
      }

      // Check if minor merger
      if ((stmax_pair.second.data().Mass / stmax_pair.first.data().Mass >= minor_merger_ratio) &&
          (stmax_pair.second.data().Mass / stmax_pair.first.data().Mass <= 1.0/minor_merger_ratio)) {
        // Check if most recent
        if (no_minor_merger_yet_baryonic) {
          md.SnapNumLastMinorMergerBaryonic[index_orig] = sub.snapnum();
          no_minor_merger_yet_baryonic = false;
        }
        // Add to merger counters
        md.NumMinorMergersTotalBaryonic[index_orig] += 1;
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 1.0)
          md.NumMinorMergersLastGyrBaryonic[index_orig] += 1;
        if (redshifts_all[sub.snapnum()] < 1.0)
          md.NumMinorMergersSinceRedshiftOneBaryonic[index_orig] += 1;
        if (redshifts_all[sub.snapnum()] < 2.0)
          md.NumMinorMergersSinceRedshiftTwoBaryonic[index_orig] += 1;
        if (redshifts_all[sub.snapnum()] < 3.0)
          md.NumMinorMergersSinceRedshiftThreeBaryonic[index_orig] += 1;
      }

      ///////////////////
      // STELLAR MASSES //
      ///////////////////

      // Only proceed if both stellar masses at stmax are > 0
      if (!((stmax_pair.first.data().SubhaloMassType[4] > 0) &&
          (stmax_pair.second.data().SubhaloMassType[4] > 0))) {
        continue;
      }

      // Check if major merger
      if ((stmax_pair.second.data().SubhaloMassType[4] / stmax_pair.first.data().SubhaloMassType[4] >= major_merger_ratio) &&
          (stmax_pair.second.data().SubhaloMassType[4] / stmax_pair.first.data().SubhaloMassType[4] <= 1.0/major_merger_ratio)) {
        // Check if most recent
        if (no_major_merger_yet) {
          md.SnapNumLastMajorMerger[index_orig] = sub.snapnum();
          no_major_merger_yet = false;
        }
        // Add to merger counters
        md.NumMajorMergersTotal[index_orig] += 1;
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 1.0)
          md.NumMajorMergersLastGyr[index_orig] += 1;
        if (redshifts_all[sub.snapnum()] < 1.0)
          md.NumMajorMergersSinceRedshiftOne[index_orig] += 1;
        if (redshifts_all[sub.snapnum()] < 2.0)
          md.NumMajorMergersSinceRedshiftTwo[index_orig] += 1;
        if (redshifts_all[sub.snapnum()] < 3.0)
          md.NumMajorMergersSinceRedshiftThree[index_orig] += 1;
      }

      // Check if minor merger
      if ((stmax_pair.second.data().SubhaloMassType[4] / stmax_pair.first.data().SubhaloMassType[4] >= minor_merger_ratio) &&
          (stmax_pair.second.data().SubhaloMassType[4] / stmax_pair.first.data().SubhaloMassType[4] <= 1.0/minor_merger_ratio)) {
        // Check if most recent
        if (no_minor_merger_yet) {
          md.SnapNumLastMinorMerger[index_orig] = sub.snapnum();
          no_minor_merger_yet = false;
        }
        // Add to merger counters
        md.NumMinorMergersTotal[index_orig] += 1;
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 1.0)
          md.NumMinorMergersLastGyr[index_orig] += 1;
        if (redshifts_all[sub.snapnum()] < 1.0)
          md.NumMinorMergersSinceRedshiftOne[index_orig] += 1;
        if (redshifts_all[sub.snapnum()] < 2.0)
          md.NumMinorMergersSinceRedshiftTwo[index_orig] += 1;
        if (redshifts_all[sub.snapnum()] < 3.0)
          md.NumMinorMergersSinceRedshiftThree[index_orig] += 1;
      }
    }

    // Next iteration
    sub = first_prog;
    first_prog = sub.first_progenitor();
  }
}

/** @brief Get snapshot of last major (minor) merger. */
void merger_history_all(const std::string& basedir, const std::string& treedir,
    const std::string& writepath, const snapnum_type snapnum_first,
    const snapnum_type snapnum_last) {

  // Get time (in Gyr) and redshift for each snapshot.
  auto redshifts_all = get_redshifts();
  auto times_all = get_times_Gyr();

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

    // Store data here:
    MergerData md(nsubs);

    // Iterate over subhalos and find the last major (minor) merger.
    auto snap = tree.snapshot(snapnum);
    for (auto sub_it = snap.begin(); sub_it != snap.end(); ++sub_it) {
      merger_history_sub(*sub_it, redshifts_all, times_all, md);
    }

    // Write to file.
    wall_clock.start();
    add_array(writefile, md.SnapNumLastMajorMerger, "SnapNumLastMajorMerger", H5::PredType::NATIVE_INT16);
    add_array(writefile, md.NumMajorMergersLastGyr, "NumMajorMergersLastGyr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMajorMergersTotal, "NumMajorMergersTotal", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMajorMergersSinceRedshiftOne, "NumMajorMergersSinceRedshiftOne", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMajorMergersSinceRedshiftTwo, "NumMajorMergersSinceRedshiftTwo", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMajorMergersSinceRedshiftThree, "NumMajorMergersSinceRedshiftThree", H5::PredType::NATIVE_UINT32);

    add_array(writefile, md.SnapNumLastMinorMerger, "SnapNumLastMinorMerger", H5::PredType::NATIVE_INT16);
    add_array(writefile, md.NumMinorMergersLastGyr, "NumMinorMergersLastGyr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMinorMergersTotal, "NumMinorMergersTotal", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMinorMergersSinceRedshiftOne, "NumMinorMergersSinceRedshiftOne", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMinorMergersSinceRedshiftTwo, "NumMinorMergersSinceRedshiftTwo", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMinorMergersSinceRedshiftThree, "NumMinorMergersSinceRedshiftThree", H5::PredType::NATIVE_UINT32);

    add_array(writefile, md.SnapNumLastMajorMergerBaryonic, "SnapNumLastMajorMergerBaryonic", H5::PredType::NATIVE_INT16);
    add_array(writefile, md.NumMajorMergersLastGyrBaryonic, "NumMajorMergersLastGyrBaryonic", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMajorMergersTotalBaryonic, "NumMajorMergersTotalBaryonic", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMajorMergersSinceRedshiftOneBaryonic, "NumMajorMergersSinceRedshiftOneBaryonic", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMajorMergersSinceRedshiftTwoBaryonic, "NumMajorMergersSinceRedshiftTwoBaryonic", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMajorMergersSinceRedshiftThreeBaryonic, "NumMajorMergersSinceRedshiftThreeBaryonic", H5::PredType::NATIVE_UINT32);

    add_array(writefile, md.SnapNumLastMinorMergerBaryonic, "SnapNumLastMinorMergerBaryonic", H5::PredType::NATIVE_INT16);
    add_array(writefile, md.NumMinorMergersLastGyrBaryonic, "NumMinorMergersLastGyrBaryonic", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMinorMergersTotalBaryonic, "NumMinorMergersTotalBaryonic", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMinorMergersSinceRedshiftOneBaryonic, "NumMinorMergersSinceRedshiftOneBaryonic", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMinorMergersSinceRedshiftTwoBaryonic, "NumMinorMergersSinceRedshiftTwoBaryonic", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMinorMergersSinceRedshiftThreeBaryonic, "NumMinorMergersSinceRedshiftThreeBaryonic", H5::PredType::NATIVE_UINT32);

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
  merger_history_all(basedir, treedir, writepath, snapnum_first, snapnum_last);

  // Print wall clock time and speedup
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  std::cout << "Speedup: " << cpu_clock.seconds()/wall_clock.seconds() << ".\n";

  return 0;
}
