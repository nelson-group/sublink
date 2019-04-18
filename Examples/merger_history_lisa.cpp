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

  // Last 2 Gyr
  std::vector<count_type> NumMajorMergersLast2Gyr;
  std::vector<count_type> NumMinorMergersLast2Gyr;
  std::vector<count_type> NumMergersLast2Gyr;
  std::vector<real_type> MeanStellarMassRatioLast2Gyr;
  std::vector<real_type> MeanStellarMassLast2Gyr;
  std::vector<real_type> MeanRedshiftAtPeakMassLast2Gyr;
  std::vector<real_type> AccretedStellarMassLast2Gyr;

  // Last 5 Gyr
  std::vector<count_type> NumMajorMergersLast5Gyr;
  std::vector<count_type> NumMinorMergersLast5Gyr;
  std::vector<count_type> NumMergersLast5Gyr;
  std::vector<real_type> MeanStellarMassRatioLast5Gyr;
  std::vector<real_type> MeanStellarMassLast5Gyr;
  std::vector<real_type> MeanRedshiftAtPeakMassLast5Gyr;
  std::vector<real_type> AccretedStellarMassLast5Gyr;

  // Last 8 Gyr
  std::vector<count_type> NumMajorMergersLast8Gyr;
  std::vector<count_type> NumMinorMergersLast8Gyr;
  std::vector<count_type> NumMergersLast8Gyr;
  std::vector<real_type> MeanStellarMassRatioLast8Gyr;
  std::vector<real_type> MeanStellarMassLast8Gyr;
  std::vector<real_type> MeanRedshiftAtPeakMassLast8Gyr;
  std::vector<real_type> AccretedStellarMassLast8Gyr;

  // Since z=5
  std::vector<count_type> NumMajorMergersSinceRedshift5;
  std::vector<count_type> NumMinorMergersSinceRedshift5;
  std::vector<count_type> NumMergersSinceRedshift5;
  std::vector<real_type> MeanStellarMassRatioSinceRedshift5;
  std::vector<real_type> MeanStellarMassSinceRedshift5;
  std::vector<real_type> MeanRedshiftAtPeakMassSinceRedshift5;
  std::vector<real_type> AccretedStellarMassSinceRedshift5;

  /** Constructor.
   * @param[in] nsubs The number of subhalos in this snapshot.
   */
  MergerData(const uint32_t nsubs)
      : NumMajorMergersLast2Gyr(std::vector<count_type>(nsubs, 0)),
        NumMinorMergersLast2Gyr(std::vector<count_type>(nsubs, 0)),
        NumMergersLast2Gyr(std::vector<count_type>(nsubs, 0)),
        MeanStellarMassRatioLast2Gyr(std::vector<real_type>(nsubs, 0)),
        MeanStellarMassLast2Gyr(std::vector<real_type>(nsubs, 0)),
        MeanRedshiftAtPeakMassLast2Gyr(std::vector<real_type>(nsubs, 0)),
        AccretedStellarMassLast2Gyr(std::vector<real_type>(nsubs, 0)),

        NumMajorMergersLast5Gyr(std::vector<count_type>(nsubs, 0)),
        NumMinorMergersLast5Gyr(std::vector<count_type>(nsubs, 0)),
        NumMergersLast5Gyr(std::vector<count_type>(nsubs, 0)),
        MeanStellarMassRatioLast5Gyr(std::vector<real_type>(nsubs, 0)),
        MeanStellarMassLast5Gyr(std::vector<real_type>(nsubs, 0)),
        MeanRedshiftAtPeakMassLast5Gyr(std::vector<real_type>(nsubs, 0)),
        AccretedStellarMassLast5Gyr(std::vector<real_type>(nsubs, 0)),

        NumMajorMergersLast8Gyr(std::vector<count_type>(nsubs, 0)),
        NumMinorMergersLast8Gyr(std::vector<count_type>(nsubs, 0)),
        NumMergersLast8Gyr(std::vector<count_type>(nsubs, 0)),
        MeanStellarMassRatioLast8Gyr(std::vector<real_type>(nsubs, 0)),
        MeanStellarMassLast8Gyr(std::vector<real_type>(nsubs, 0)),
        MeanRedshiftAtPeakMassLast8Gyr(std::vector<real_type>(nsubs, 0)),
        AccretedStellarMassLast8Gyr(std::vector<real_type>(nsubs, 0)),

        NumMajorMergersSinceRedshift5(std::vector<count_type>(nsubs, 0)),
        NumMinorMergersSinceRedshift5(std::vector<count_type>(nsubs, 0)),
        NumMergersSinceRedshift5(std::vector<count_type>(nsubs, 0)),
        MeanStellarMassRatioSinceRedshift5(std::vector<real_type>(nsubs, 0)),
        MeanStellarMassSinceRedshift5(std::vector<real_type>(nsubs, 0)),
        MeanRedshiftAtPeakMassSinceRedshift5(std::vector<real_type>(nsubs, 0)),
        AccretedStellarMassSinceRedshift5(std::vector<real_type>(nsubs, 0)) {
  }
};

/** @brief Calculate merging history for a given subhalo. */
void merger_history_sub(Subhalo sub,
    const std::vector<real_type>& redshifts_all,
    const std::vector<real_type>& times_all,
    MergerData& md) {

  uint32_t index_orig = sub.index();
  snapnum_type snapnum_orig = sub.snapnum();

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

      // Check if major merger
      if (mass_ratio >= major_merger_ratio) {
        // Add to merger counters
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 2.0)
          md.NumMajorMergersLast2Gyr[index_orig] += 1;
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 5.0)
          md.NumMajorMergersLast5Gyr[index_orig] += 1;
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 8.0)
          md.NumMajorMergersLast8Gyr[index_orig] += 1;
        if (redshifts_all[sub.snapnum()] < 5.0)
          md.NumMajorMergersSinceRedshift5[index_orig] += 1;
      }

      // Check if minor merger (does not overlap with major mergers)
      if ((mass_ratio >= minor_merger_ratio) &&
           (mass_ratio <  major_merger_ratio)) {
        // Add to merger counters
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 2.0)
          md.NumMinorMergersLast2Gyr[index_orig] += 1;
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 5.0)
          md.NumMinorMergersLast5Gyr[index_orig] += 1;
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 8.0)
          md.NumMinorMergersLast8Gyr[index_orig] += 1;
        if (redshifts_all[sub.snapnum()] < 5.0)
          md.NumMinorMergersSinceRedshift5[index_orig] += 1;
      }

      // All mergers
      if (true) {
        // Add to merger counters
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 2.0)
          md.NumMergersLast2Gyr[index_orig] += 1;
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 5.0)
          md.NumMergersLast5Gyr[index_orig] += 1;
        if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 8.0)
          md.NumMergersLast8Gyr[index_orig] += 1;
        if (redshifts_all[sub.snapnum()] < 5.0)
          md.NumMergersSinceRedshift5[index_orig] += 1;
      }

      // Add contributions
      // Last 2 Gyr
      if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 2.0) {
        md.MeanStellarMassRatioLast2Gyr[index_orig] += smallest_mstar * mass_ratio;
        md.MeanStellarMassLast2Gyr[index_orig] += smallest_mstar * smallest_mstar;
        md.MeanRedshiftAtPeakMassLast2Gyr[index_orig] += smallest_mstar * redshifts_all[snapnum_stmax];
        md.AccretedStellarMassLast2Gyr[index_orig] += smallest_mstar;
      }
      // Last 5 Gyr
      if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 5.0) {
        md.MeanStellarMassRatioLast5Gyr[index_orig] += smallest_mstar * mass_ratio;
        md.MeanStellarMassLast5Gyr[index_orig] += smallest_mstar * smallest_mstar;
        md.MeanRedshiftAtPeakMassLast5Gyr[index_orig] += smallest_mstar * redshifts_all[snapnum_stmax];
        md.AccretedStellarMassLast5Gyr[index_orig] += smallest_mstar;
      }
      // Last 8 Gyr
      if (times_all[snapnum_orig] - times_all[sub.snapnum()] < 8.0) {
        md.MeanStellarMassRatioLast8Gyr[index_orig] += smallest_mstar * mass_ratio;
        md.MeanStellarMassLast8Gyr[index_orig] += smallest_mstar * smallest_mstar;
        md.MeanRedshiftAtPeakMassLast8Gyr[index_orig] += smallest_mstar * redshifts_all[snapnum_stmax];
        md.AccretedStellarMassLast8Gyr[index_orig] += smallest_mstar;
      }
      // Since z=5
      if (redshifts_all[sub.snapnum()] < 5.0) {
        md.MeanStellarMassRatioSinceRedshift5[index_orig] += smallest_mstar * mass_ratio;
        md.MeanStellarMassSinceRedshift5[index_orig] += smallest_mstar * smallest_mstar;
        md.MeanRedshiftAtPeakMassSinceRedshift5[index_orig] += smallest_mstar * redshifts_all[snapnum_stmax];
        md.AccretedStellarMassSinceRedshift5[index_orig] += smallest_mstar;
      }
    }

    // Next iteration
    sub = first_prog;
    first_prog = sub.first_progenitor();
  }

  // Divide by denominator to get an average
  // Last 2 Gyr
  if (md.AccretedStellarMassLast2Gyr[index_orig] > 0) {
    real_type denominator = md.AccretedStellarMassLast2Gyr[index_orig];
    md.MeanStellarMassRatioLast2Gyr[index_orig] /= denominator;
    md.MeanStellarMassLast2Gyr[index_orig] /= denominator;
    md.MeanRedshiftAtPeakMassLast2Gyr[index_orig] /= denominator;
  }
  else {  // quantities are undefined
    md.MeanStellarMassRatioLast2Gyr[index_orig] = -1;
    md.MeanStellarMassLast2Gyr[index_orig] = -1;
    md.MeanRedshiftAtPeakMassLast2Gyr[index_orig] = -1;
  }
  // Last 5 Gyr
  if (md.AccretedStellarMassLast5Gyr[index_orig] > 0) {
    real_type denominator = md.AccretedStellarMassLast5Gyr[index_orig];
    md.MeanStellarMassRatioLast5Gyr[index_orig] /= denominator;
    md.MeanStellarMassLast5Gyr[index_orig] /= denominator;
    md.MeanRedshiftAtPeakMassLast5Gyr[index_orig] /= denominator;
  }
  else {  // quantities are undefined
    md.MeanStellarMassRatioLast5Gyr[index_orig] = -1;
    md.MeanStellarMassLast5Gyr[index_orig] = -1;
    md.MeanRedshiftAtPeakMassLast5Gyr[index_orig] = -1;
  }
  // Last 8 Gyr
  if (md.AccretedStellarMassLast8Gyr[index_orig] > 0) {
    real_type denominator = md.AccretedStellarMassLast8Gyr[index_orig];
    md.MeanStellarMassRatioLast8Gyr[index_orig] /= denominator;
    md.MeanStellarMassLast8Gyr[index_orig] /= denominator;
    md.MeanRedshiftAtPeakMassLast8Gyr[index_orig] /= denominator;
  }
  else {  // quantities are undefined
    md.MeanStellarMassRatioLast8Gyr[index_orig] = -1;
    md.MeanStellarMassLast8Gyr[index_orig] = -1;
    md.MeanRedshiftAtPeakMassLast8Gyr[index_orig] = -1;
  }
  // Since z=5
  if (md.AccretedStellarMassSinceRedshift5[index_orig] > 0) {
    real_type denominator = md.AccretedStellarMassSinceRedshift5[index_orig];
    md.MeanStellarMassRatioSinceRedshift5[index_orig] /= denominator;
    md.MeanStellarMassSinceRedshift5[index_orig] /= denominator;
    md.MeanRedshiftAtPeakMassSinceRedshift5[index_orig] /= denominator;
  }
  else {  // quantities are undefined
    md.MeanStellarMassRatioSinceRedshift5[index_orig] = -1;
    md.MeanStellarMassSinceRedshift5[index_orig] = -1;
    md.MeanRedshiftAtPeakMassSinceRedshift5[index_orig] = -1;
  }
}

/** @brief Get snapshot of last major (minor) merger. */
void merger_history_all(
    const std::string& suite, const std::string& basedir,
    const std::string& treedir, const std::string& writepath,
    const snapnum_type snapnum_first, const snapnum_type snapnum_last) {

  (void) suite;  // Silence unused variable warning

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

    // Iterate over subhalos and do merger history calculations.
    auto snap = tree.snapshot(snapnum);
    for (auto sub_it = snap.begin(); sub_it != snap.end(); ++sub_it) {
      merger_history_sub(*sub_it, redshifts_all, times_all, md);
    }

    // Write to file.
    wall_clock.start();

    add_array(writefile, md.NumMajorMergersLast2Gyr, "NumMajorMergersLast2Gyr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMinorMergersLast2Gyr, "NumMinorMergersLast2Gyr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMergersLast2Gyr, "NumMergersLast2Gyr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.MeanStellarMassRatioLast2Gyr, "MeanStellarMassRatioLast2Gyr", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, md.MeanStellarMassLast2Gyr, "MeanStellarMassLast2Gyr", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, md.MeanRedshiftAtPeakMassLast2Gyr, "MeanRedshiftAtPeakMassLast2Gyr", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, md.AccretedStellarMassLast2Gyr, "AccretedStellarMassLast2Gyr", H5::PredType::NATIVE_FLOAT);

    add_array(writefile, md.NumMajorMergersLast5Gyr, "NumMajorMergersLast5Gyr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMinorMergersLast5Gyr, "NumMinorMergersLast5Gyr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMergersLast5Gyr, "NumMergersLast5Gyr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.MeanStellarMassRatioLast5Gyr, "MeanStellarMassRatioLast5Gyr", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, md.MeanStellarMassLast5Gyr, "MeanStellarMassLast5Gyr", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, md.MeanRedshiftAtPeakMassLast5Gyr, "MeanRedshiftAtPeakMassLast5Gyr", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, md.AccretedStellarMassLast5Gyr, "AccretedStellarMassLast5Gyr", H5::PredType::NATIVE_FLOAT);

    add_array(writefile, md.NumMajorMergersLast8Gyr, "NumMajorMergersLast8Gyr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMinorMergersLast8Gyr, "NumMinorMergersLast8Gyr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMergersLast8Gyr, "NumMergersLast8Gyr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.MeanStellarMassRatioLast8Gyr, "MeanStellarMassRatioLast8Gyr", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, md.MeanStellarMassLast8Gyr, "MeanStellarMassLast8Gyr", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, md.MeanRedshiftAtPeakMassLast8Gyr, "MeanRedshiftAtPeakMassLast8Gyr", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, md.AccretedStellarMassLast8Gyr, "AccretedStellarMassLast8Gyr", H5::PredType::NATIVE_FLOAT);

    add_array(writefile, md.NumMajorMergersSinceRedshift5, "NumMajorMergersSinceRedshift5", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMinorMergersSinceRedshift5, "NumMinorMergersSinceRedshift5", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.NumMergersSinceRedshift5, "NumMergersSinceRedshift5", H5::PredType::NATIVE_UINT32);
    add_array(writefile, md.MeanStellarMassRatioSinceRedshift5, "MeanStellarMassRatioSinceRedshift5", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, md.MeanStellarMassSinceRedshift5, "MeanStellarMassSinceRedshift5", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, md.MeanRedshiftAtPeakMassSinceRedshift5, "MeanRedshiftAtPeakMassSinceRedshift5", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, md.AccretedStellarMassSinceRedshift5, "AccretedStellarMassSinceRedshift5", H5::PredType::NATIVE_FLOAT);

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
