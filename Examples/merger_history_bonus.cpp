/** @file merger_history_bonus.cpp
 * @brief Calculate a few more merger statistics, such as the
 *        mean (mass-weighted) mass ratio of all mergers, the mean
 *        (mass-weighted) time since all mergers, and the mean
 *        (mass-weighted) cold gas fraction of the secondary progenitors.
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
#include "../Util/Cosmology.hpp"  // cosmo::get_redshifts, cosmo::get_times_Gyr
#include "../Util/TreeUtil.hpp"

/** @brief Type for number of mergers. */
typedef uint32_t count_type;

struct MergerData {
  // The following quantities are weighted by the stellar mass
  // of the secondary progenitor:
  std::vector<real_type> MeanMassRatio;
  std::vector<real_type> MeanRedshift;
  std::vector<real_type> MeanLookbackTime;
  std::vector<real_type> MeanGasFraction;  // star-forming gas

  /** Constructor.
   * @param[in] nsubs The number of subhalos in this snapshot.
   */
  MergerData(const uint32_t nsubs)
      : MeanMassRatio(std::vector<real_type>(nsubs, -1)),
        MeanRedshift(std::vector<real_type>(nsubs, -1)),
        MeanLookbackTime(std::vector<real_type>(nsubs, -1)),
        MeanGasFraction(std::vector<real_type>(nsubs, -1)) {
  }
};

/** @brief Find the last major/minor mergers for a given subhalo. */
void merger_history_sub(Subhalo sub,
    const std::vector<real_type>& redshifts_all,
    const std::vector<real_type>& times_all,
    MergerData& md) {

  auto sub_orig = sub;

  // Sum of the stellar masses of all secondary progenitors:
  real_type denominator = 0;

  // Initialize mean quantities to zero:
  md.MeanMassRatio[sub_orig.index()] = 0;
  md.MeanRedshift[sub_orig.index()] = 0;
  md.MeanLookbackTime[sub_orig.index()] = 0;
  md.MeanGasFraction[sub_orig.index()] = 0;

  // Iterate over first progenitor
  auto first_prog = sub_orig.first_progenitor();
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
      real_type mstar_stmax_1 = stmax_pair.first.data().SubhaloMassType[4];
      real_type mstar_stmax_2 = stmax_pair.second.data().SubhaloMassType[4];
      if (!((mstar_stmax_1 > 0) && (mstar_stmax_2 > 0))) {
        continue;
      }

      // Make sure that mass ratio < 1:
      real_type mass_ratio = std::min(mstar_stmax_2/mstar_stmax_1,
          mstar_stmax_1/mstar_stmax_2);

      // "Weights" also correspond to the smallest of the two masses:
      real_type smallest_mstar = std::min(mstar_stmax_1, mstar_stmax_2);

      // Gas fraction is an "intensive" quantity: should always take the secondary
      real_type mgal_stmax_2 = stmax_pair.second.data().Mass;
      real_type fgas_stmax_2 = (mgal_stmax_2 - mstar_stmax_2) / mgal_stmax_2;

      // Add contributions
      md.MeanMassRatio[sub_orig.index()] += smallest_mstar * mass_ratio;
      md.MeanRedshift[sub_orig.index()] += smallest_mstar * redshifts_all[sub.snapnum()];
      md.MeanLookbackTime[sub_orig.index()] += smallest_mstar * (times_all[sub_orig.snapnum()] - times_all[sub.snapnum()]);
      md.MeanGasFraction[sub_orig.index()] += smallest_mstar * fgas_stmax_2;
      denominator += smallest_mstar;
    }

    // Next iteration
    sub = first_prog;
    first_prog = sub.first_progenitor();
  }

  // Divide by denominator to get an average:
  if (denominator > 0) {
    md.MeanMassRatio[sub_orig.index()] /= denominator;
    md.MeanRedshift[sub_orig.index()] /= denominator;
    md.MeanLookbackTime[sub_orig.index()] /= denominator;
    md.MeanGasFraction[sub_orig.index()] /= denominator;
  }
  else {  // quantities are undefined
    md.MeanMassRatio[sub_orig.index()] = -1;
    md.MeanRedshift[sub_orig.index()] = -1;
    md.MeanLookbackTime[sub_orig.index()] = -1;
    md.MeanGasFraction[sub_orig.index()] = -1;

  }
}

/** @brief Get snapshot of last major (minor) merger. */
void merger_history_all(const std::string& basedir, const std::string& treedir,
    const std::string& writepath, const snapnum_type snapnum_first,
    const snapnum_type snapnum_last) {

  // Get time (in Gyr) and redshift for each snapshot.
  auto redshifts_all = cosmo::get_redshifts();
  auto times_all = cosmo::get_times_Gyr();

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

    // Iterate over subhalos and calculate some merger statistics.
    auto snap = tree.snapshot(snapnum);
    for (auto sub_it = snap.begin(); sub_it != snap.end(); ++sub_it) {
      merger_history_sub(*sub_it, redshifts_all, times_all, md);
    }

    // Write to file.
    add_array(writefile, md.MeanMassRatio, "MeanMassRatio", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, md.MeanRedshift, "MeanRedshift", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, md.MeanLookbackTime, "MeanLookbackTime", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, md.MeanGasFraction, "MeanGasFraction", H5::PredType::NATIVE_FLOAT);
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
