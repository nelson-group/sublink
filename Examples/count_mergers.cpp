/** @file count_mergers.cpp
 * @brief Count mergers and print their main properties into files.
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
#include "../Util/TreeUtil.hpp"

/** @brief Datatype to store merger info. */
struct MergerInfo{
  // Descendant properties.
  real_type mstar_0;
  real_type m200_0;
  real_type delta_0;
  int8_t is_central_0;

  // Progenitor properties right before merger.
  real_type mstar_instant_1;
  real_type mstar_instant_2;

  // Progenitor properties at stellar tmax.
  real_type mstar_stmax_1;
  real_type mstar_stmax_2;
  snapnum_type snapnum_stmax;

  // Progenitor properties at infall.
  real_type mstar_infall_1;
  real_type mstar_infall_2;
  snapnum_type snapnum_infall;

  /** Constructor. */
  MergerInfo(
      real_type mstar_0_,
      real_type m200_0_,
      real_type delta_0_,
      int8_t is_central_0_,
      real_type mstar_instant_1_,
      real_type mstar_instant_2_,
      real_type mstar_stmax_1_,
      real_type mstar_stmax_2_,
      snapnum_type snapnum_stmax_,
      real_type mstar_infall_1_,
      real_type mstar_infall_2_,
      snapnum_type snapnum_infall_)
      : mstar_0(mstar_0_),
        m200_0(m200_0_),
        delta_0(delta_0_),
        is_central_0(is_central_0_),
        mstar_instant_1(mstar_instant_1_),
        mstar_instant_2(mstar_instant_2_),
        mstar_stmax_1(mstar_stmax_1_),
        mstar_stmax_2(mstar_stmax_2_),
        snapnum_stmax(snapnum_stmax_),
        mstar_infall_1(mstar_infall_1_),
        mstar_infall_2(mstar_infall_2_),
        snapnum_infall(snapnum_infall_) {
  }

  /** Disable default constructor. */
  MergerInfo() = delete;
};

/** @brief Corresponding HDF5 datatype. */
static H5::CompType H5MergerInfo() {
  H5::CompType mtype(sizeof(MergerInfo));
  // Descendant properties
  mtype.insertMember("mstar_0", HOFFSET(MergerInfo, mstar_0), H5::PredType::NATIVE_FLOAT);
  mtype.insertMember("m200_0", HOFFSET(MergerInfo, m200_0), H5::PredType::NATIVE_FLOAT);
  mtype.insertMember("delta_0", HOFFSET(MergerInfo, delta_0), H5::PredType::NATIVE_FLOAT);
  // Although 1 byte is enough, this entry might occupy 8 bytes:
  mtype.insertMember("is_central_0", HOFFSET(MergerInfo, is_central_0), H5::PredType::NATIVE_INT8);
  // Progenitor properties right before merger
  mtype.insertMember("mstar_instant_1", HOFFSET(MergerInfo, mstar_instant_1), H5::PredType::NATIVE_FLOAT);
  mtype.insertMember("mstar_instant_2", HOFFSET(MergerInfo, mstar_instant_2), H5::PredType::NATIVE_FLOAT);
  // Progenitor properties at stellar tmax
  mtype.insertMember("mstar_stmax_1", HOFFSET(MergerInfo, mstar_stmax_1), H5::PredType::NATIVE_FLOAT);
  mtype.insertMember("mstar_stmax_2", HOFFSET(MergerInfo, mstar_stmax_2), H5::PredType::NATIVE_FLOAT);
  mtype.insertMember("snapnum_stmax", HOFFSET(MergerInfo, snapnum_stmax), H5::PredType::NATIVE_INT16);
  // Progenitor properties at infall
  mtype.insertMember("mstar_infall_1", HOFFSET(MergerInfo, mstar_infall_1), H5::PredType::NATIVE_FLOAT);
  mtype.insertMember("mstar_infall_2", HOFFSET(MergerInfo, mstar_infall_2), H5::PredType::NATIVE_FLOAT);
  mtype.insertMember("snapnum_infall", HOFFSET(MergerInfo, snapnum_infall), H5::PredType::NATIVE_INT16);
  return mtype;
}


/** @brief Count mergers for a given subhalo and print data to file.
 * @param[in] sub The subhalo of interest.
 * @param[in] writefile An open file to write the results.
 */
void count_mergers_sub(Subhalo sub, std::vector<MergerInfo>& merger_info,
    const std::vector<std::vector<real_type>>& overdensities) {
  // Only proceed if @a sub has at least one progenitor.
  auto first_prog = sub.first_progenitor();
  if (!first_prog.is_valid())
    return;

  // Iterate over "next progenitor" link to count mergers.
  for (auto next_prog = first_prog.next_progenitor();
      next_prog.is_valid(); next_prog = next_prog.next_progenitor()) {

    // Descendant properties.
    real_type mstar_0 = sub.data().SubhaloMassType[4];
    real_type m200_0 = sub.data().Group_M_Crit200;
    real_type delta_0 = overdensities[sub.data().SnapNum][sub.data().SubfindID];
    auto is_central_0 = static_cast<int8_t>(sub.data().FirstSubhaloInFOFGroupID == sub.data().SubhaloID);

    // Progenitor properties right before merger.
    real_type mstar_instant_1 = first_prog.data().SubhaloMassType[4];
    real_type mstar_instant_2 = next_prog.data().SubhaloMassType[4];

    // Progenitor properties at stellar tmax.
    real_type mstar_stmax_1 = -1;
    real_type mstar_stmax_2 = -1;
    snapnum_type snapnum_stmax = -1;
    auto stmax_pair = get_stmax_pair(first_prog, next_prog);
    if (stmax_pair.first.is_valid()) {
      mstar_stmax_1 = stmax_pair.first.data().SubhaloMassType[4];
      mstar_stmax_2 = stmax_pair.second.data().SubhaloMassType[4];
      snapnum_stmax = stmax_pair.second.data().SnapNum;
    }

    // Progenitor properties at infall.
    real_type mstar_infall_1 = -1;
    real_type mstar_infall_2 = -1;
    snapnum_type snapnum_infall = -1;
    auto infall_pair = get_infall_pair(first_prog, next_prog);
    if (infall_pair.first.is_valid()) {
      mstar_infall_1 = infall_pair.first.data().SubhaloMassType[4];
      mstar_infall_2 = infall_pair.second.data().SubhaloMassType[4];
      snapnum_infall = infall_pair.second.data().SnapNum;
    }

    // Add to data structure
    auto cur_merger_info = MergerInfo(mstar_0, m200_0, delta_0, is_central_0,
        mstar_instant_1, mstar_instant_2, mstar_stmax_1, mstar_stmax_2, snapnum_stmax,
        mstar_infall_1, mstar_infall_2, snapnum_infall);

    merger_info.push_back(cur_merger_info);
  }
}

/** @brief Count mergers and print to files. */
void count_mergers_all(const std::string& simdir,
    const std::string& envdir, const std::string& treedir,
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

    // Environment filename
    std::stringstream tmp_stream;
    tmp_stream << envdir << "/environment_" <<
        std::setfill('0') << std::setw(3) << snapnum << ".hdf5";
    std::string filename = tmp_stream.str();

    // Only proceed if file exists
    std::ifstream file(filename);
    if (!file) {
      // No overdensities to read. Fill array with -999's
      overdensities[snapnum].resize(nsubs, -999);
      continue;
    }
    file.close();

    // If file is empty, fill array with -999's
    auto h5file = H5::H5File(filename, H5F_ACC_RDONLY);
    if (!H5Lexists(h5file.getId(), "delta", H5P_DEFAULT)) {
      h5file.close();
      overdensities[snapnum].resize(nsubs, -999);
      continue;
    }
    h5file.close();

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

    // Store merger info here
    std::vector<MergerInfo> merger_info;

    // Iterate over subhalos and count mergers.
    for (auto sub_it = snap.begin(); sub_it != snap.end(); ++sub_it) {
      count_mergers_sub(*sub_it, merger_info, overdensities);
    }

    // Filename for this snapshot.
    std::stringstream tmp_stream;
    tmp_stream << writepath << "_" <<
        std::setfill('0') << std::setw(3) << snapnum << ".hdf5";
    std::string writefilename = tmp_stream.str();

    // Write to HDF5 file.
    H5::H5File writefile(writefilename, H5F_ACC_TRUNC);
    add_array(writefile, merger_info, "Mergers", H5MergerInfo());
    writefile.close();

    std::cout << "Finished for snapshot " << snapnum << std::endl;
  }
}


int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 7) {
    std::cerr << "Usage: " << argv[0] << " simdir envdir treedir writepath" <<
        " snapnum_first snapnum_last\n";
    exit(1);
  }

  // Read input
  std::string simdir(argv[1]);
  std::string envdir(argv[2]);
  std::string treedir(argv[3]);
  std::string writepath(argv[4]);
  snapnum_type snapnum_first = atoi(argv[5]);
  snapnum_type snapnum_last = atoi(argv[6]);

  // Measure CPU and wall clock (real) time
  WallClock wall_clock;
  CPUClock cpu_clock;

  // Do stuff
  count_mergers_all(simdir, envdir, treedir, writepath, snapnum_first, snapnum_last);

  // Print wall clock time and speedup
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";
  std::cout << "Speedup: " << cpu_clock.seconds()/wall_clock.seconds() << ".\n";

  return 0;
}
