/** @file find_descendants.cpp
 * @brief Find subhalo descendants for a given range of snapshots.
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

#include <string>

#include "ParticleMatcher.hpp"


/** Get some types from ParticleMatcher and make them our own. */
typedef typename ParticleMatcher::snapnum_type snapnum_type;

int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 10) {
    std::cerr << "Usage: " << argv[0] << " basedir writepath " <<
        "snapnum_first snapnum_last snapnum_start snapnum_end " <<
        "tracking_scheme pass skipsnaps_filename";
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

  // Measure CPU and wall clock (real) time
  CPUClock cpu_clock;
  WallClock wall_clock;

  // Do stuff
  auto pm = ParticleMatcher(basedir, basedir, snapnum_last-1, snapnum_last,
      tracking_scheme);

  // Print CPU and wall clock time
  std::cout << "Finished.\n";
  std::cout << "CPU time: "  << cpu_clock.seconds() << " s.\n";
  std::cout << "Wall clock time: "  << wall_clock.seconds() << " s.\n";

  return 0;
}
