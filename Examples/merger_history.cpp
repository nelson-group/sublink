/** @file merger_history.cpp
 * @brief For each galaxy, print some simple merger statistics, such as
 *        the snapshot number of the last major merger and the number of
 *        major mergers in the last Gyr.
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
  std::vector<snapnum_type> SnapNumNextMajorMerger;
  std::vector<count_type> NumMajorMergersLast250Myr;
  std::vector<count_type> NumMajorMergersLast500Myr;
  std::vector<count_type> NumMajorMergersLastGyr;
  std::vector<count_type> NumMajorMergersSinceRedshiftOne;
  std::vector<count_type> NumMajorMergersSinceRedshiftTwo;
  std::vector<count_type> NumMajorMergersTotal;

  // Minor mergers
  std::vector<snapnum_type> SnapNumLastMinorMerger;
  std::vector<snapnum_type> SnapNumNextMinorMerger;
  std::vector<count_type> NumMinorMergersLast250Myr;
  std::vector<count_type> NumMinorMergersLast500Myr;
  std::vector<count_type> NumMinorMergersLastGyr;
  std::vector<count_type> NumMinorMergersSinceRedshiftOne;
  std::vector<count_type> NumMinorMergersSinceRedshiftTwo;
  std::vector<count_type> NumMinorMergersTotal;

  // All mergers
  std::vector<snapnum_type> SnapNumLastMerger;
  std::vector<snapnum_type> SnapNumNextMerger;
  std::vector<count_type> NumMergersLast250Myr;
  std::vector<count_type> NumMergersLast500Myr;
  std::vector<count_type> NumMergersLastGyr;
  std::vector<count_type> NumMergersSinceRedshiftOne;
  std::vector<count_type> NumMergersSinceRedshiftTwo;
  std::vector<count_type> NumMergersTotal;

  /** Constructor.
   * @param[in] nsubs The number of subhalos in this snapshot.
   */
  MergerData(const uint32_t nsubs)
      : SnapNumLastMajorMerger(std::vector<snapnum_type>(nsubs, -1)),
        SnapNumNextMajorMerger(std::vector<snapnum_type>(nsubs, -1)),
        NumMajorMergersLast250Myr(std::vector<count_type>(nsubs, 0)),
        NumMajorMergersLast500Myr(std::vector<count_type>(nsubs, 0)),
        NumMajorMergersLastGyr(std::vector<count_type>(nsubs, 0)),
        NumMajorMergersSinceRedshiftOne(std::vector<count_type>(nsubs, 0)),
        NumMajorMergersSinceRedshiftTwo(std::vector<count_type>(nsubs, 0)),
        NumMajorMergersTotal(std::vector<count_type>(nsubs, 0)),

        SnapNumLastMinorMerger(std::vector<snapnum_type>(nsubs, -1)),
        SnapNumNextMinorMerger(std::vector<snapnum_type>(nsubs, -1)),
        NumMinorMergersLast250Myr(std::vector<count_type>(nsubs, 0)),
        NumMinorMergersLast500Myr(std::vector<count_type>(nsubs, 0)),
        NumMinorMergersLastGyr(std::vector<count_type>(nsubs, 0)),
        NumMinorMergersSinceRedshiftOne(std::vector<count_type>(nsubs, 0)),
        NumMinorMergersSinceRedshiftTwo(std::vector<count_type>(nsubs, 0)),
        NumMinorMergersTotal(std::vector<count_type>(nsubs, 0)),

        SnapNumLastMerger(std::vector<snapnum_type>(nsubs, -1)),
        SnapNumNextMerger(std::vector<snapnum_type>(nsubs, -1)),
        NumMergersLast250Myr(std::vector<count_type>(nsubs, 0)),
        NumMergersLast500Myr(std::vector<count_type>(nsubs, 0)),
        NumMergersLastGyr(std::vector<count_type>(nsubs, 0)),
        NumMergersSinceRedshiftOne(std::vector<count_type>(nsubs, 0)),
        NumMergersSinceRedshiftTwo(std::vector<count_type>(nsubs, 0)),
        NumMergersTotal(std::vector<count_type>(nsubs, 0)) {
  }
};

/** @brief Find the last major/minor mergers for a given subhalo. */
void merger_history_sub(Subhalo sub,
    const std::vector<real_type>& redshifts_all,
    const std::vector<real_type>& times_all,
    MergerData& md) {

  auto sub_orig = sub;
  uint32_t index_orig = sub.index();
  snapnum_type snapnum_orig = sub.snapnum();

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
          no_major_merger_yet = false;
        }
        // Add to merger counters
        md.NumMajorMergersTotal[index_orig] += 1;
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 0.25) {
          md.NumMajorMergersLast250Myr[index_orig] += 1;
        }
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 0.5) {
          md.NumMajorMergersLast500Myr[index_orig] += 1;
        }
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 1.0) {
          md.NumMajorMergersLastGyr[index_orig] += 1;
        }
        if (redshifts_all[sub.snapnum()] < 1.0) {
          md.NumMajorMergersSinceRedshiftOne[index_orig] += 1;
        }
        if (redshifts_all[sub.snapnum()] < 2.0) {
          md.NumMajorMergersSinceRedshiftTwo[index_orig] += 1;
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
          no_minor_merger_yet = false;
        }
        // Add to merger counters
        md.NumMinorMergersTotal[index_orig] += 1;
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 0.25) {
          md.NumMinorMergersLast250Myr[index_orig] += 1;
        }
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 0.5) {
          md.NumMinorMergersLast500Myr[index_orig] += 1;
        }
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 1.0) {
          md.NumMinorMergersLastGyr[index_orig] += 1;
        }
        if (redshifts_all[sub.snapnum()] < 1.0) {
          md.NumMinorMergersSinceRedshiftOne[index_orig] += 1;
        }
        if (redshifts_all[sub.snapnum()] < 2.0) {
          md.NumMinorMergersSinceRedshiftTwo[index_orig] += 1;
        }
      }

      // All mergers
      if (true) {
        // Check if most recent
        if (no_merger_yet) {
          md.SnapNumLastMerger[index_orig] = sub.snapnum();
          no_merger_yet = false;
        }
        // Add to merger counters
        md.NumMergersTotal[index_orig] += 1;
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 0.25) {
          md.NumMergersLast250Myr[index_orig] += 1;
        }
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 0.5) {
          md.NumMergersLast500Myr[index_orig] += 1;
        }
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 1.0) {
          md.NumMergersLastGyr[index_orig] += 1;
        }
        if (redshifts_all[sub.snapnum()] < 1.0) {
          md.NumMergersSinceRedshiftOne[index_orig] += 1;
        }
        if (redshifts_all[sub.snapnum()] < 2.0) {
          md.NumMergersSinceRedshiftTwo[index_orig] += 1;
        }
      }
    }

    // Next iteration
    sub = first_prog;
    first_prog = sub.first_progenitor();
  }

  // Now iterate forward in time (to find time until next merger)
  no_major_merger_yet = true;
  no_minor_merger_yet = true;
  no_merger_yet = true;

  // Get descendant of original subhalo
  sub = sub_orig.descendant();
  while (sub.is_valid()) {
    // Get first progenitor of descendant
    first_prog = sub.first_progenitor();

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
        // Check if earliest
        if (no_major_merger_yet) {
          md.SnapNumNextMajorMerger[index_orig] = sub.snapnum();
          no_major_merger_yet = false;
        }
      }

      // Check if minor merger (does not overlap with major mergers)
      if (((mass_ratio >= minor_merger_ratio) &&
           (mass_ratio <  major_merger_ratio)) ||
          ((mass_ratio >  1.0/major_merger_ratio) &&
           (mass_ratio <= 1.0/minor_merger_ratio))) {
        // Check if earliest
        if (no_minor_merger_yet) {
          md.SnapNumNextMinorMerger[index_orig] = sub.snapnum();
          no_minor_merger_yet = false;
        }
      }

      // All mergers
      if (true) {
        // Check if earliest
        if (no_merger_yet) {
          md.SnapNumNextMerger[index_orig] = sub.snapnum();
          no_merger_yet = false;
        }
      }
    }

    // Next iteration
    sub = sub.descendant();
  }

}

/** @brief Get snapshot of last major (minor) merger. */
void merger_history_all(
    const std::string& suite, const std::string& basedir,
    const std::string& treedir, const std::string& writepath,
    const snapnum_type snapnum_first, const snapnum_type snapnum_last) {

  // Get time (in Gyr) and redshift for each snapshot.
  auto redshifts_all = cosmo::get_redshifts(suite);
  auto times_all = cosmo::get_times_Gyr(suite);

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
      merger_history_sub(*sub_it, redshifts_all, times_all, md);
    }

    // Write to file.
    wall_clock.start();

    add_array(writefile, md.SnapNumLastMajorMerger, "SnapNumLastMajorMerger", H5::PredType::NATIVE_INT16);
    add_array(writefile, md.SnapNumNextMajorMerger, "SnapNumNextMajorMerger", H5::PredType::NATIVE_INT16);
    add_array(writefile, md.NumMajorMergersLast250Myr, "NumMajorMergersLast250Myr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMajorMergersLast500Myr, "NumMajorMergersLast500Myr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMajorMergersLastGyr, "NumMajorMergersLastGyr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMajorMergersSinceRedshiftOne, "NumMajorMergersSinceRedshiftOne", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMajorMergersSinceRedshiftTwo, "NumMajorMergersSinceRedshiftTwo", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMajorMergersTotal, "NumMajorMergersTotal", H5::PredType::NATIVE_UINT32);

    add_array(writefile, md.SnapNumLastMinorMerger, "SnapNumLastMinorMerger", H5::PredType::NATIVE_INT16);
    add_array(writefile, md.SnapNumNextMinorMerger, "SnapNumNextMinorMerger", H5::PredType::NATIVE_INT16);
    add_array(writefile, md.NumMinorMergersLast250Myr, "NumMinorMergersLast250Myr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMinorMergersLast500Myr, "NumMinorMergersLast500Myr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMinorMergersLastGyr, "NumMinorMergersLastGyr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMinorMergersSinceRedshiftOne, "NumMinorMergersSinceRedshiftOne", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMinorMergersSinceRedshiftTwo, "NumMinorMergersSinceRedshiftTwo", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMinorMergersTotal, "NumMinorMergersTotal", H5::PredType::NATIVE_UINT32);

    add_array(writefile, md.SnapNumLastMerger, "SnapNumLastMerger", H5::PredType::NATIVE_INT16);
    add_array(writefile, md.SnapNumNextMerger, "SnapNumNextMerger", H5::PredType::NATIVE_INT16);
    add_array(writefile, md.NumMergersLast250Myr, "NumMergersLast250Myr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMergersLast500Myr, "NumMergersLast500Myr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMergersLastGyr, "NumMergersLastGyr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMergersSinceRedshiftOne, "NumMergersSinceRedshiftOne", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMergersSinceRedshiftTwo, "NumMergersSinceRedshiftTwo", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMergersTotal, "NumMergersTotal", H5::PredType::NATIVE_UINT32);

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

  // Print wall clock time and speedup
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  std::cout << "Speedup: " << cpu_clock.seconds()/wall_clock.seconds() << ".\n";

  return 0;
}
