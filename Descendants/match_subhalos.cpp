/** @file match_subhalos.cpp
 * @brief Find subhalo matches between two different simulations.
 *
 * @author Vicente Rodriguez-Gomez (v.rodriguez@irya.unam.mx)
 */

#include <fstream>

#include "ParticleMatcher.hpp"

int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 10) {
    std::cerr << "Usage: " << argv[0] << " basedir1 basedir2 writepath " <<
        "snapnum_first snapnum_last snapnum_start snapnum_end " <<
        "tracking_scheme skipsnaps_filename alpha_weight\n";
    exit(1);
  }

  // Read input
  std::string basedir1(argv[1]);
  std::string basedir2(argv[2]);
  std::string writepath(argv[3]);
  snapnum_type snapnum_first = atoi(argv[4]);
  snapnum_type snapnum_last = atoi(argv[5]);
  snapnum_type snapnum_start = atoi(argv[6]);
  snapnum_type snapnum_end = atoi(argv[7]);
  std::string tracking_scheme(argv[8]);  // "Subhalos" or "Galaxies"
  std::string skipsnaps_filename(argv[9]);
  real_type alpha_weight = atof(argv[10]);  // Usually 0 or 1

  // Create list of valid snapshots
  auto valid_snapnums = get_valid_snapnums(skipsnaps_filename,
      snapnum_first, snapnum_last);

  // Iterate over snapshot range
  for (auto snapnum1 = snapnum_start; snapnum1 <= snapnum_end; ++snapnum1) {
    // Check that first snapshot is valid
    auto it = std::find(valid_snapnums.begin(), valid_snapnums.end(), snapnum1);
    if (it == valid_snapnums.end())
      continue;

    // Define second snapshot number
    snapnum_type snapnum2 = snapnum1;

    // Measure CPU and wall clock (real) time
    CPUClock cpu_clock;
    WallClock wall_clock;

    // Find descendants and write to files
    auto pm = ParticleMatcher(basedir1, basedir2, snapnum1, snapnum2,
                              tracking_scheme, alpha_weight);
    pm.write_to_file(writepath, false);

    // Print CPU and wall clock time
    std::cout << "Finished.\n";
    std::cout << "CPU time: "  << cpu_clock.seconds() << " s.\n";
    std::cout << "Wall clock time: "  << wall_clock.seconds() << " s.\n";
    std::cout << "\n";
  }

  return 0;
}
