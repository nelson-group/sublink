/** @file find_descendants.cpp
 * @brief Find subhalo descendants for a given range of snapshots.
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

#include <fstream>

#include "ParticleMatcher.hpp"

int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 10) {
    std::cerr << "Usage: " << argv[0] << " basedir writepath " <<
        "snapnum_first snapnum_last snapnum_start snapnum_end " <<
        "tracking_scheme pass skipsnaps_filename\n";
    exit(1);
  }

  // Read input
  std::string basedir(argv[1]);
  std::string writepath(argv[2]);
  snapnum_type snapnum_first = atoi(argv[3]);
  snapnum_type snapnum_last = atoi(argv[4]);
  snapnum_type snapnum_start = atoi(argv[5]);
  snapnum_type snapnum_end = atoi(argv[6]);
  std::string tracking_scheme(argv[7]);  /* Subhalos or Galaxies */
  std::string pass(argv[8]);  /* first or second  */
  std::string skipsnaps_filename(argv[9]);

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
    snapnum_type snapnum2 = -1;
    if (pass == "first") {
      if (it+1 < valid_snapnums.end())
        snapnum2 = *(it+1);
    }
    else if (pass == "second") {
      if (it+2 < valid_snapnums.end())
        snapnum2 = *(it+2);
    }
    else
      assert(false);

    // Measure CPU and wall clock (real) time
    CPUClock cpu_clock;
    WallClock wall_clock;

    // Find descendants
    auto pm = ParticleMatcher(basedir, basedir, snapnum1, snapnum2, tracking_scheme);
    pm.write_to_file(writepath);

    // Print CPU and wall clock time
    std::cout << "Finished.\n";
    std::cout << "CPU time: "  << cpu_clock.seconds() << " s.\n";
    std::cout << "Wall clock time: "  << wall_clock.seconds() << " s.\n";
    std::cout << "\n";
  }

  return 0;
}
