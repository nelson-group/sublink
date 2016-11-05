/** @file build_trees.cpp
 * @brief Construct merger trees using descendant files.
 *
 * @author Vicente Rodriguez-Gomez (vrg@jhu.edu)
 */

#include "TreesClasses.hpp"

int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 6) {
    std::cerr << "Usage: ./BuildTrees input_path output_path snapnum_first " <<
        "snapnum_last skipsnaps_filename\n";
    exit(1);
  }

  // Read input
  std::string input_path(argv[1]);
  std::string output_path(argv[2]);
  snapnum_type snapnum_first = atoi(argv[3]);
  snapnum_type snapnum_last = atoi(argv[4]);
  std::string skipsnaps_filename(argv[5]);

  // Measure CPU and wall clock (real) time
  CPUClock cpu_clock;
  WallClock wall_clock;

  // Construct trees and write to files.
  auto all_trees = AllTrees(input_path, output_path, snapnum_first,
      snapnum_last, skipsnaps_filename);
  (void) all_trees;  // silence compiler warning

  // Print CPU and wall clock time
  std::cout << "Finished.\n";
  std::cout << "CPU time: "  << cpu_clock.seconds() << " s.\n";
  std::cout << "Wall clock time: "  << wall_clock.seconds() << " s.\n";
  std::cout << "\n";

  return 0;
}
