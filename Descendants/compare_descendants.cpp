/** @file compare_descendants.cpp
 * @brief Compare the descendants created by @a find_descendants.cpp
 * during the first and second "passes."
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

#include <fstream>

#include "../InputOutput/GeneralHDF5.hpp"
#include "ParticleMatcher.hpp"

/** @brief Compare the descendants from the first and second "passes."
 *
 * For each subhalo at snapshot1, we compare the "skipped" descendant
 * at snapshot3, obtained by skipping snapshot2, with the "straight" descendant
 * at snapshot3, a.k.a. the "descendant of the descendant". If the two possible
 * descendants at snapshot3 are the same object, we define the final one as the
 * "immediate" one at snapshot2. Otherwise, we define the unique descendant
 * as the "skipped" one at snapshot3, since it is the one with the highest
 * score.
 *
 * Update (01/23/15): Also keep track of the second highest score
 * at each snapshot.
 */
void compare_descendants(const snapnum_type snapnum1,
    const snapnum_type snapnum2, const std::string& writepath,
    const bool trivial) {

  // Create filenames
  std::stringstream tmp_stream;
  tmp_stream << writepath << "_first_" <<
                std::setfill('0') << std::setw(3) << snapnum1 << ".hdf5";
  std::string filename_12 = tmp_stream.str();
  tmp_stream.str("");
  tmp_stream << writepath << "_second_" <<
                std::setfill('0') << std::setw(3) << snapnum1 << ".hdf5";
  std::string filename_13 = tmp_stream.str();
  tmp_stream.str("");
  tmp_stream << writepath << "_first_" <<
                std::setfill('0') << std::setw(3) << snapnum2 << ".hdf5";
  std::string filename_23 = tmp_stream.str();
  tmp_stream.str("");
  tmp_stream << writepath << "_" <<
                std::setfill('0') << std::setw(3) << snapnum1 << ".hdf5";
  std::string writefilename = tmp_stream.str();

  // Open output file
  H5::H5File writefile(writefilename, H5F_ACC_TRUNC);

  // Only proceed if first file is non-empty
  H5::H5File tmp_file(filename_12, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (!H5Lexists(tmp_file.getId(), "/DescendantIndex", H5P_DEFAULT)) {
    std::cout << "Snapshot " << snapnum1 << ": Skipping empty file.\n";
    writefile.close();
    tmp_file.close();
    return;
  }
  tmp_file.close();

  // Read info from first snapshot
  auto sub_len_12 = read_dataset<uint32_t>(filename_12, "SubhaloLen");
  auto sub_mass_12 = read_dataset<real_type>(filename_12, "SubhaloMass");
  auto sub_grnr_12 = read_dataset<uint32_t>(filename_12, "SubhaloGrNr");
  auto desc_index_12 = read_dataset<index_type>(filename_12, "DescendantIndex");
  auto first_score_12 = read_dataset<real_type>(filename_12, "FirstScore");
  auto second_score_12 = read_dataset<real_type>(filename_12, "SecondScore");
  uint32_t nsubs = sub_len_12.size();

  // Indicate if snapshot 2 is skipped
  std::vector<uint8_t> skip_snapshot(nsubs, 0);

  // If we cannot skip snapshots, not much to do here (trivial case)
  if (trivial) {
    add_array(writefile, sub_len_12, "SubhaloLen", H5::PredType::NATIVE_UINT32);
    add_array(writefile, sub_mass_12, "SubhaloMass", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, sub_grnr_12, "SubhaloGrNr", H5::PredType::NATIVE_UINT32);
    add_array(writefile, desc_index_12, "DescendantIndex", H5::PredType::NATIVE_INT32);
    add_array(writefile, first_score_12, "FirstScore", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, second_score_12, "SecondScore", H5::PredType::NATIVE_FLOAT);
    add_array(writefile, skip_snapshot, "SkipSnapshot", H5::PredType::NATIVE_UINT8);
    writefile.close();
    return;
  }

  // Check that the other descendant files are non-empty
  H5::H5File tmp_file_13(filename_13, H5F_ACC_RDONLY, H5P_DEFAULT);
  H5::H5File tmp_file_23(filename_23, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (!H5Lexists(tmp_file_13.getId(), "/DescendantIndex", H5P_DEFAULT) ||
      !H5Lexists(tmp_file_23.getId(), "/DescendantIndex", H5P_DEFAULT)) {
    std::cerr << "BAD: Missing some descendant files.\n";
    writefile.close();
    tmp_file_13.close();
    tmp_file_23.close();
    return;
  }
  tmp_file_13.close();
  tmp_file_23.close();

  // Read info from other files
  auto desc_index_13 = read_dataset<index_type>(filename_13, "DescendantIndex");
  auto desc_index_23 = read_dataset<index_type>(filename_23, "DescendantIndex");
  auto first_score_13 = read_dataset<real_type>(filename_13, "FirstScore");
  auto second_score_13 = read_dataset<real_type>(filename_13, "SecondScore");

  // Initialize descendants and scores to the ones from snapshot 2
  auto desc_final = desc_index_12;
  auto first_score_final = first_score_12;
  auto second_score_final = second_score_12;

  // Compare descendants
  for (uint32_t sub_index = 0; sub_index < nsubs; ++sub_index) {
    auto desc_immediate = desc_index_12[sub_index];
    auto desc_skip = desc_index_13[sub_index];
    auto first_score_skip = first_score_13[sub_index];
    auto second_score_skip = second_score_13[sub_index];

    if (desc_skip == -1)
      continue;

    if (desc_immediate == -1) {
      desc_final[sub_index] = desc_skip;
      first_score_final[sub_index] = first_score_skip;
      second_score_final[sub_index] = second_score_skip;
      skip_snapshot[sub_index] = 1;
      continue;
    }

    auto desc_straight = desc_index_23[desc_immediate];
    if (desc_straight != desc_skip) {
      desc_final[sub_index] = desc_skip;
      first_score_final[sub_index] = first_score_skip;
      second_score_final[sub_index] = second_score_skip;
      skip_snapshot[sub_index] = 1;
    }
  }

  // Write to file
  add_array(writefile, sub_len_12, "SubhaloLen", H5::PredType::NATIVE_UINT32);
  add_array(writefile, sub_mass_12, "SubhaloMass", H5::PredType::NATIVE_FLOAT);
  add_array(writefile, sub_grnr_12, "SubhaloGrNr", H5::PredType::NATIVE_UINT32);
  add_array(writefile, desc_final, "DescendantIndex", H5::PredType::NATIVE_INT32);
  add_array(writefile, first_score_final, "FirstScore", H5::PredType::NATIVE_FLOAT);
  add_array(writefile, second_score_final, "SecondScore", H5::PredType::NATIVE_FLOAT);
  add_array(writefile, skip_snapshot, "SkipSnapshot", H5::PredType::NATIVE_UINT8);
  writefile.close();
}


int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0] << " writepath " <<
        "snapnum_first snapnum_last skipsnaps_filename\n";
    exit(1);
  }

  // Read input
  std::string writepath(argv[1]);
  snapnum_type snapnum_first = atoi(argv[2]);
  snapnum_type snapnum_last = atoi(argv[3]);
  std::string skipsnaps_filename(argv[4]);

  // Create list of valid snapshots
  auto valid_snapnums = get_valid_snapnums(skipsnaps_filename,
      snapnum_first, snapnum_last);

  // Iterate over snapshot range
  for (auto snapnum1 = snapnum_first; snapnum1 <= snapnum_last; ++snapnum1) {
    // Check that first snapshot is valid
    auto it = std::find(valid_snapnums.begin(), valid_snapnums.end(), snapnum1);
    if (it == valid_snapnums.end())
      continue;

    // Define next two snapshots
    snapnum_type snapnum2 = -1;
    snapnum_type snapnum3 = -1;
    if (it+1 < valid_snapnums.end())
      snapnum2 = *(it+1);
    if (it+2 < valid_snapnums.end())
      snapnum3 = *(it+2);

    // Use "trivial" label in cases when we cannot skip snapshots
    bool trivial = false;
    if ((snapnum2 == -1) || (snapnum3 == -1))
      trivial = true;

    // Measure wall clock (real) time
    WallClock wall_clock;

    // Compare descendants
    compare_descendants(snapnum1, snapnum2, writepath, trivial);

    // Print wall clock time
    std::cout << "Finished for snapshot " << snapnum1 << ".\n";
    std::cout << "Time: "  << wall_clock.seconds() << " s.\n";
    std::cout << "\n";
  }

  return 0;
}
