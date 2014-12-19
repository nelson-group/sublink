#include <iostream>
#include <chrono>

#include "../InputOutput/ReadArepoHDF5.hpp"

int main(int argc, char** argv)
{
  // Check input arguments
  if (argc != 6) {
    std::cerr << "Usage: " << argv[0] << " basedir treedir writepath" <<
        " snapnum_first snapnum_last\n";
    exit(1);
  }

//  // Read input
//  std::string basedir(argv[1]);
//  std::string treedir(argv[2]);
//  std::string writepath(argv[3]);
//  int16_t snapnum_first = atoi(argv[4]);
//  int16_t snapnum_last = atoi(argv[5]);

  // Measure real time
  auto start = std::chrono::system_clock::now();

  // Test: read attribute from header
  std::string filename("/n/ghernquist/Illustris/Runs/L75n455FP/output/snapdir_135/snap_135.0.hdf5");

  auto box_size = arepo_get_attribute<double>(filename, "BoxSize");

  std::cout << "Box size is " << box_size << " kpc/h.\n";

  // Print elapsed time
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Elapsed time: " << elapsed_seconds.count() << " seconds.\n";

  return 0;
}
