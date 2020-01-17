/** @file infall_catalog.cpp
 * @brief For each halo (FoF group), get information about all the objects
 *        that have ever "infalled" into it.
 *
 * @author Vicente Rodriguez-Gomez (v.rodriguez@irya.unam.mx)
 */

// Include some extra quantities from the merger trees:
#define INFALL_CATALOG

#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "../InputOutput/ReadArepoHDF5.hpp"
#include "../InputOutput/ReadSubfindHDF5.hpp"
#include "../InputOutput/ReadTreeHDF5.hpp"
#include "../InputOutput/GeneralHDF5.hpp"
#include "../Util/GeneralUtil.hpp"
#include "../Util/TreeUtil.hpp"

/** @brief Get some info for all the objects that have "infalled"
 * into FoF groups at z=0.
 */
void infall_catalog(const std::string& basedir, const std::string& treedir,
    const std::string& writedir, const snapnum_type snapnum_last) {

  // Get box size in ckpc/h; note type conversion
  std::stringstream tmp_stream;
  tmp_stream << basedir << "/snapdir_" <<
      std::setfill('0') << std::setw(3) << snapnum_last << "/snap_" <<
      std::setfill('0') << std::setw(3) << snapnum_last << ".0.hdf5";
  std::string file_name = tmp_stream.str();
  real_type box_size = static_cast<real_type>(
      arepo::get_scalar_attribute<double>(file_name, "BoxSize"));

  // Load merger tree
  std::cout << "Loading merger tree...\n";
  WallClock wall_clock;
  int filenum = -1;  // read from all merger tree files
  std::string name = "tree_extended";  // Full format
  Tree tree(treedir, name, filenum);
  std::cout << "Loaded merger tree. Total time: " <<
      wall_clock.seconds() << " s.\n";
  std::cout << "\n";

  // Load some FoF group info at z=0
  std::cout << "Loading FoF group info...\n";
  wall_clock.start();
  auto group_first_sub = subfind::read_block<uint32_t>(
      basedir, snapnum_last, "Group", "GroupFirstSub", -1);
  auto ngroups = subfind::get_scalar_attribute<uint32_t>(
      basedir, snapnum_last, "Ngroups_Total");
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";

  // Store info on-the-fly into these arrays
  std::vector<index_type> GroupIndex;
  // Properties at infall:
  std::vector<snapnum_type> SnapNumInfall;
  std::vector<index_type> SubfindIDInfall;
  std::vector<real_type> SubhaloMassInfall;
  std::vector<FloatArray<6>> SubhaloMassTypeInfall;
  std::vector<real_type> SubhaloVmaxInfall;
  // Properties at the latest appearance of each subhalo:
  std::vector<snapnum_type> SnapNumLastIdentified;
  std::vector<index_type> SubfindIDLastIdentified;
  std::vector<real_type> SubhaloMassLastIdentified;
  std::vector<FloatArray<6>> SubhaloMassTypeLastIdentified;
  std::vector<real_type> SubhaloVmaxLastIdentified;
  // Properties at last virial crossing (R200):
  std::vector<snapnum_type> SnapNumLastCrossing;
  std::vector<index_type> SubfindIDLastCrossing;
  std::vector<real_type> SubhaloMassLastCrossing;
  std::vector<FloatArray<6>> SubhaloMassTypeLastCrossing;
  std::vector<real_type> SubhaloVmaxLastCrossing;
  // Properties at first virial crossing (R200):
  std::vector<snapnum_type> SnapNumFirstCrossing;
  std::vector<index_type> SubfindIDFirstCrossing;
  std::vector<real_type> SubhaloMassFirstCrossing;
  std::vector<FloatArray<6>> SubhaloMassTypeFirstCrossing;
  std::vector<real_type> SubhaloVmaxFirstCrossing;

  // Store offsets here
  std::vector<uint32_t> GroupFirstSub(ngroups, 0);
  std::vector<uint32_t> GroupNsubs(ngroups, 0);

  // For cases where infall is well defined but R200Crit is not:
  FloatArray<6> invalid_masstype = {-1, -1, -1, -1, -1, -1};

  // Iterate over FoF groups at z=0
  std::cout << "Iterating over FoF groups...\n";
  wall_clock.start();
  for (uint32_t group_index = 0; group_index < ngroups; ++group_index) {

    // Offset of this FoF group
    if (group_index > 0)
      GroupFirstSub[group_index] = GroupFirstSub[group_index-1] + GroupNsubs[group_index-1];

    // Count subhalos that "infalled" onto this FoF group.
    uint32_t subhalo_count = 0;

    // If this FoF group has no subhalos, skip.
    auto first_sub_index = group_first_sub[group_index];
    if (first_sub_index == static_cast<uint32_t>(-1))
      continue;

    // Iterate over first progenitor of central, which we call
    // orig_sub at every snapshot.
    for (auto orig_sub = tree.subhalo(snapnum_last, first_sub_index);
        orig_sub.is_valid(); orig_sub = orig_sub.first_progenitor()) {

      // Iterate over subhalos in same FoF group as orig_sub.
      for (auto cur_sub = orig_sub.first_subhalo_in_fof_group();
          cur_sub.is_valid(); cur_sub = cur_sub.next_subhalo_in_fof_group()) {

        // orig_sub is not always the first subhalo in the FoF
        if (cur_sub == orig_sub)
          continue;

        // If cur_sub is the first progenitor of its descendant, we
        // already took it into account, so we skip.
        auto desc = cur_sub.descendant();
        if (desc.is_valid() && desc.first_progenitor() == cur_sub)
          continue;

        // Get properties at infall, taking orig_sub as the primary
        auto infall_pair = get_infall_pair(orig_sub, cur_sub);

        // Only proceed if infall is well defined
        auto snapnum_infall = infall_pair.second.snapnum();
        if (snapnum_infall == -1)
          continue;

        // If got here, we have found a new, valid "infall object".
        GroupIndex.push_back(group_index);

        SnapNumInfall.push_back(snapnum_infall);
        SubfindIDInfall.push_back(infall_pair.second.index());
        SubhaloMassInfall.push_back(infall_pair.second.data().SubhaloMass);
        SubhaloMassTypeInfall.push_back(infall_pair.second.data().SubhaloMassType);
        SubhaloVmaxInfall.push_back(infall_pair.second.data().SubhaloVmax);

        SnapNumLastIdentified.push_back(cur_sub.snapnum());
        SubfindIDLastIdentified.push_back(cur_sub.index());
        SubhaloMassLastIdentified.push_back(cur_sub.data().SubhaloMass);
        SubhaloMassTypeLastIdentified.push_back(cur_sub.data().SubhaloMassType);
        SubhaloVmaxLastIdentified.push_back(cur_sub.data().SubhaloVmax);

        // Add info about moment of last (latest) virial crossing
        auto pair_last_crossing = get_pair_virial_crossing(orig_sub, cur_sub,
            box_size, true);
        auto snapnum_last_crossing = pair_last_crossing.second.snapnum();
        if (snapnum_last_crossing == -1) {
          SnapNumLastCrossing.push_back(-1);
          SubfindIDLastCrossing.push_back(-1);
          SubhaloMassLastCrossing.push_back(-1);
          SubhaloMassTypeLastCrossing.push_back(invalid_masstype);
          SubhaloVmaxLastCrossing.push_back(-1);
        }
        else {
          SnapNumLastCrossing.push_back(snapnum_last_crossing);
          SubfindIDLastCrossing.push_back(pair_last_crossing.second.index());
          SubhaloMassLastCrossing.push_back(pair_last_crossing.second.data().SubhaloMass);
          SubhaloMassTypeLastCrossing.push_back(pair_last_crossing.second.data().SubhaloMassType);
          SubhaloVmaxLastCrossing.push_back(pair_last_crossing.second.data().SubhaloVmax);
        }

        // Add info about moment of first (earliest) virial crossing
        auto pair_first_crossing = get_pair_virial_crossing(orig_sub, cur_sub,
            box_size, false);
        auto snapnum_first_crossing = pair_first_crossing.second.snapnum();
        if (snapnum_first_crossing == -1) {
          SnapNumFirstCrossing.push_back(-1);
          SubfindIDFirstCrossing.push_back(-1);
          SubhaloMassFirstCrossing.push_back(-1);
          SubhaloMassTypeFirstCrossing.push_back(invalid_masstype);
          SubhaloVmaxFirstCrossing.push_back(-1);
        }
        else {
          SnapNumFirstCrossing.push_back(snapnum_first_crossing);
          SubfindIDFirstCrossing.push_back(pair_first_crossing.second.index());
          SubhaloMassFirstCrossing.push_back(pair_first_crossing.second.data().SubhaloMass);
          SubhaloMassTypeFirstCrossing.push_back(pair_first_crossing.second.data().SubhaloMassType);
          SubhaloVmaxFirstCrossing.push_back(pair_first_crossing.second.data().SubhaloVmax);
        }

        // Add to counter
        ++subhalo_count;
      }
    }

    // Number of subhalos that "infalled" onto this FoF group
    GroupNsubs[group_index] = subhalo_count;

    // Print some output occasionally
    if (group_index+1 % 1000 == 0) {
      std::cout << "Already processed " << group_index+1 << "FoF groups.\n";
    }
  }
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";

  // Output filename
  tmp_stream.str("");
  tmp_stream << writedir << "/infall_catalog_" <<
      std::setfill('0') << std::setw(3) << snapnum_last << ".hdf5";
  std::string writefilename = tmp_stream.str();

  // Write to output file
  std::cout << "Writing to HDF5 file...\n";
  wall_clock.start();
  H5::H5File writefile(writefilename, H5F_ACC_TRUNC);
  add_array(writefile, GroupIndex, "GroupIndex", H5::PredType::NATIVE_INT32);

  add_array(writefile, SnapNumInfall, "SnapNumInfall", H5::PredType::NATIVE_INT16);
  add_array(writefile, SubfindIDInfall, "SubfindIDInfall", H5::PredType::NATIVE_INT32);
  add_array(writefile, SubhaloMassInfall, "SubhaloMassInfall", H5::PredType::NATIVE_FLOAT);
  add_array_2d(writefile, SubhaloMassTypeInfall, "SubhaloMassTypeInfall", H5::PredType::NATIVE_FLOAT);
  add_array(writefile, SubhaloVmaxInfall, "SubhaloVmaxInfall", H5::PredType::NATIVE_FLOAT);

  add_array(writefile, SnapNumLastIdentified, "SnapNumLastIdentified", H5::PredType::NATIVE_INT16);
  add_array(writefile, SubfindIDLastIdentified, "SubfindIDLastIdentified", H5::PredType::NATIVE_INT32);
  add_array(writefile, SubhaloMassLastIdentified, "SubhaloMassLastIdentified", H5::PredType::NATIVE_FLOAT);
  add_array_2d(writefile, SubhaloMassTypeLastIdentified, "SubhaloMassTypeLastIdentified", H5::PredType::NATIVE_FLOAT);
  add_array(writefile, SubhaloVmaxLastIdentified, "SubhaloVmaxLastIdentified", H5::PredType::NATIVE_FLOAT);

  add_array(writefile, SnapNumLastCrossing, "SnapNumLastCrossing", H5::PredType::NATIVE_INT16);
  add_array(writefile, SubfindIDLastCrossing, "SubfindIDLastCrossing", H5::PredType::NATIVE_INT32);
  add_array(writefile, SubhaloMassLastCrossing, "SubhaloMassLastCrossing", H5::PredType::NATIVE_FLOAT);
  add_array_2d(writefile, SubhaloMassTypeLastCrossing, "SubhaloMassTypeLastCrossing", H5::PredType::NATIVE_FLOAT);
  add_array(writefile, SubhaloVmaxLastCrossing, "SubhaloVmaxLastCrossing", H5::PredType::NATIVE_FLOAT);

  add_array(writefile, SnapNumFirstCrossing, "SnapNumFirstCrossing", H5::PredType::NATIVE_INT16);
  add_array(writefile, SubfindIDFirstCrossing, "SubfindIDFirstCrossing", H5::PredType::NATIVE_INT32);
  add_array(writefile, SubhaloMassFirstCrossing, "SubhaloMassFirstCrossing", H5::PredType::NATIVE_FLOAT);
  add_array_2d(writefile, SubhaloMassTypeFirstCrossing, "SubhaloMassTypeFirstCrossing", H5::PredType::NATIVE_FLOAT);
  add_array(writefile, SubhaloVmaxFirstCrossing, "SubhaloVmaxFirstCrossing", H5::PredType::NATIVE_FLOAT);

  writefile.close();
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";

  // Offsets filename
  tmp_stream.str("");
  tmp_stream << writedir << "/offsets_" <<
      std::setfill('0') << std::setw(3) << snapnum_last << ".hdf5";
  std::string offsets_filename = tmp_stream.str();

  // Write to offsets file
  std::cout << "Writing offsets...\n";
  wall_clock.start();
  H5::H5File offsets_file(offsets_filename, H5F_ACC_TRUNC);
  add_array(offsets_file, GroupFirstSub, "GroupFirstSub", H5::PredType::NATIVE_UINT32);
  add_array(offsets_file, GroupNsubs, "GroupNsubs", H5::PredType::NATIVE_UINT32);
  offsets_file.close();
  std::cout << "Time: " << wall_clock.seconds() << " s.\n";
}

int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0] << " basedir treedir writedir" <<
        " snapnum_last\n";
    exit(1);
  }

  // Read input
  std::string basedir(argv[1]);
  std::string treedir(argv[2]);
  std::string writedir(argv[3]);
  snapnum_type snapnum_last = atoi(argv[4]);

  // Measure CPU and wall clock (real) time
  WallClock wall_clock;
  CPUClock cpu_clock;

  // Do stuff
  infall_catalog(basedir, treedir, writedir, snapnum_last);

  // Print wall clock time and speedup
  std::cout << "Total time: " << wall_clock.seconds() << " s.\n";
  std::cout << "Speedup: " << cpu_clock.seconds()/wall_clock.seconds() << ".\n";

  return 0;
}
