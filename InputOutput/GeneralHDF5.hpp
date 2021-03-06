#pragma once

/** @file GeneralHDF5.hpp
 * @brief General-purpose utilities to read from and write to HDF5 files.
 *
 * @author Vicente Rodriguez-Gomez (v.rodriguez@irya.unam.mx)
 */

#include <vector>
#include <string>
#include <cassert>

#include "H5Cpp_wrapper.hpp"

/** @brief Check whether an HDF5 file exists or not.
 * @param[in] file_name Path to the input file.
 * @return True if HDF5 file exists; false otherwise.
 */
bool h5_file_exists(const std::string& file_name) {
    H5::Exception::dontPrint();
    try {
        H5::H5File file(file_name, H5F_ACC_RDONLY );
        file.close();
        return true;
    } catch(H5::FileIException &file_exists_err) {
        return false;
    }
    assert(false); // If this point is reached, something went wrong.
    return false;  // Return statement just to remove IDE warning.
}

/** @brief Function to read a dataset (a.k.a. block) from a single HDF5 file.
 *
 * @tparam T Type of the elements in the dataset.
 * @param[in] file_name Path to the input file.
 * @param[in] block_name Name of the dataset.
 * @return A vector with the dataset values.
 *
 * @note The return value is a vector (with an appropriately chosen type).
 *       When reading from a dataset with ndims > 1, the size() of the
 *       output vector equals the size of the first dimension of the dataset.
 */
template <typename T>
std::vector<T> read_dataset(const std::string& file_name,
    const std::string& block_name) {

  // Open dataset
  H5::H5File file(file_name, H5F_ACC_RDONLY );
  H5::DataSet dataset(file.openDataSet(block_name));

  // Get dimensions of the dataset
  H5::DataSpace file_space = dataset.getSpace();
  const unsigned int file_rank = file_space.getSimpleExtentNdims();
  hsize_t file_dims[file_rank];
  file_space.getSimpleExtentDims(file_dims, NULL);

  // Check that datatype sizes match (necessary but not sufficient)
  auto dt = dataset.getDataType();
  std::size_t row_size = dt.getSize();
  for (unsigned int k = 1; k < file_rank; ++k)
    row_size = row_size * file_dims[k];
  if (row_size != sizeof(T)) {
    std::cerr << "ERROR: mismatched datatype sizes (" << row_size <<
        " vs " << sizeof(T) << ") when reading " << block_name << ".\n";
    assert(false);
  }

  // Dataspace in memory is the same as in the HDF5 file.
  const unsigned int mem_rank = file_rank;
  hsize_t mem_dims[mem_rank];
  for (hsize_t k = 0; k < file_rank; ++k) mem_dims[k] = file_dims[k];
  H5::DataSpace mem_space(mem_rank, mem_dims);

  // Read data
  std::vector<T> retval(mem_dims[0]);  // Output vector is 1D.
  dataset.read(retval.data(), dataset.getDataType(), mem_space, file_space);

  file.close();
  return retval;
}

/** @brief Function to read a dataset (a.k.a. block) from a sequence of
 * HDF5 files with filenames ending in .0.hdf5, .1.hdf5, etc.
 *
 * @tparam T Type of the elements in the dataset.
 * @param[in] file_base Path to the input files, excluding ".N.hdf5".
 * @param[in] block_name Name of the dataset.
 * @return A vector with the dataset values.
 *
 * @note The return value is a vector (with an appropriately chosen type).
 *       When reading from a dataset with ndims > 1, the size() of the
 *       output vector equals the size of the first dimension of the
 *       (concatenated) dataset.
 */
template <typename T>
std::vector<T> read_dataset_manyfiles(const std::string& file_base,
    const std::string& block_name) {

  // Iterate over files to get number of subhalos
  std::vector<uint64_t> file_nsubs;
  std::vector<uint64_t> file_offsets;
  int64_t nsubs_in_trees = 0;
  unsigned int filenum = 0; // cannot have filenum == -1 in this function
  while (true) {
    // Path of current merger tree file
    std::stringstream tmp_stream;
    tmp_stream << file_base << "." << filenum << ".hdf5";
    std::string file_name = tmp_stream.str();

    // Only proceed if file exists
    if (!h5_file_exists(file_name))
      break;

    // Open dataset
    H5::H5File file(file_name, H5F_ACC_RDONLY );
    H5::DataSet dataset(file.openDataSet(block_name));

    // Get dimensions of the dataset (in "file space")
    H5::DataSpace file_space = dataset.getSpace();
    const unsigned int file_rank = file_space.getSimpleExtentNdims();
    hsize_t file_dims[file_rank];
    file_space.getSimpleExtentDims(file_dims, NULL);

    // Store number of subhalos
    file_nsubs.push_back(file_dims[0]);
    nsubs_in_trees += file_dims[0];
    if (filenum > 0)
      file_offsets.push_back(file_offsets[filenum-1] + file_nsubs[filenum-1]);
    else
      file_offsets.push_back(0);

    // Next iteration
    file.close();
    filenum += 1;
  }
  assert(file_nsubs.size() == file_offsets.size());
  const unsigned int nfiles = file_nsubs.size();

  // Store output in 1D vector
  std::vector<T> retval(nsubs_in_trees);

  // Iterate over files to read the data
  std::cout << "Reading " << block_name << " from " << nfiles << " files...\n";
  for (filenum = 0; filenum < nfiles; ++filenum) {
    // Path of current merger tree file
    std::stringstream tmp_stream;
    tmp_stream << file_base << "." << filenum << ".hdf5";
    std::string file_name = tmp_stream.str();

    // Open dataset (the file must exist)
    H5::H5File file(file_name, H5F_ACC_RDONLY);
    H5::DataSet dataset(file.openDataSet(block_name));

    // Get dimensions of the dataset (in "file space")
    H5::DataSpace file_space = dataset.getSpace();
    const unsigned int file_rank = file_space.getSimpleExtentNdims();
    hsize_t file_dims[file_rank];
    file_space.getSimpleExtentDims(file_dims, NULL);
    assert(file_dims[0] == file_nsubs[filenum]); // sanity check

    // Check that datatype sizes match (necessary but not sufficient)
    auto dt = dataset.getDataType();
    std::size_t row_size = dt.getSize();
    for (unsigned int k = 1; k < file_rank; ++k)
      row_size = row_size * file_dims[k];
    if (row_size != sizeof(T)) {
      std::cerr << "ERROR: mismatched datatype sizes (" << row_size <<
          " vs " << sizeof(T) << ") when reading " << block_name << ".\n";
      assert(false);
    }

    // Define the "memory dataspace", which corresponds to the full,
    // concatenated data. This is the same for all files.
    const unsigned int mem_rank = file_rank;
    hsize_t mem_dims[mem_rank];
    for (hsize_t k = 0; k < file_rank; ++k)
      mem_dims[k] = file_dims[k];
    mem_dims[0] = nsubs_in_trees; // note this line in particular
    H5::DataSpace mem_space(mem_rank, mem_dims);

    // Define "memory hyperslab" (see readdata.cpp from HDF5 tutorial)
    hsize_t offset_out[mem_rank]; // hyperslab offset
    hsize_t count_out[mem_rank]; // hyperslab size
    for (hsize_t k = 0; k < mem_rank; ++k) {
      offset_out[k] = 0;
      count_out[k] = mem_dims[k];
    }
    offset_out[0] = file_offsets[filenum];
    count_out[0] = file_nsubs[filenum];
    mem_space.selectHyperslab(H5S_SELECT_SET, count_out, offset_out);

    // Read data
    dataset.read(retval.data(), dataset.getDataType(), mem_space, file_space);

    // Close HDF5 file
    file.close();
  }
  return retval;
}

/** @brief Convenience function. Calls @a read_dataset when @a filenum >= 0
 * or @a read_dataset_manyfiles when @a filenum == -1.
 *
 * @tparam T Type of the elements in the dataset.
 * @param[in] file_base Path to the input files, excluding ".N.hdf5".
 * @param[in] filenum File number (-1 to read sequence of files).
 * @param[in] block_name Name of the dataset.
 * @return A vector with the dataset values.
 */
template <typename T>
std::vector<T> read_dataset_by_filenum(const std::string& file_base,
    const int filenum, const std::string& block_name) {
  // If given a valid file number, only read data from that file.
  if (filenum >= 0) {
    std::stringstream tmp_stream;
    tmp_stream << file_base << "." << filenum << ".hdf5";
    std::string file_name = tmp_stream.str();
    return read_dataset<T>(file_name, block_name);
  }
  // Otherwise read data from all files.
  if (filenum != -1) {
    std::cerr << "WARNING: you should use filenum == -1 when reading " <<
        "from many files.\n";
  }
  return read_dataset_manyfiles<T>(file_base, block_name);
}

/** @brief Function to add a new one-dimensional array to an open HDF5 file.
 *
 * @note Byte order is little-endian by default.
 */
template <typename T>
void add_array(H5::H5File& file, const std::vector<T>& array,
    const std::string& array_name, H5::DataType datatype) {

  // Only proceed if array is non-empty
  if (array.size() == 0)
    return;

  // Define (one-dimensional) dataspace
  hsize_t dimsf[1];  // dataset dimensions
  dimsf[0] = array.size();
  H5::DataSpace dataspace(1, dimsf);  // rank == 1

  // Create dataset
  H5::DataSet dataset = file.createDataSet(array_name, datatype, dataspace);

  // Write to dataset using default memory space
  dataset.write(array.data(), datatype, dataspace);
}

/** @brief Function to add a new two-dimensional array to an open HDF5 file.
 * @tparam T Must be the FloatArray type defined in TreeTypes.hpp,
 *           or at least something that defines a size().
 * @note Byte order is little-endian by default.
 */
template <typename T>
void add_array_2d(H5::H5File& file, const std::vector<T>& array,
    const std::string& array_name, H5::DataType datatype) {

  // Only proceed if array is non-empty
  if (array.size() == 0)
    return;

  // Define (two-dimensional) dataspace
  hsize_t dimsf[2];  // dataset dimensions
  dimsf[0] = array.size();
  dimsf[1] = sizeof(T) / sizeof(typename T::value_type);

  H5::DataSpace dataspace(2, dimsf);  // rank == 2

  // Create dataset
  H5::DataSet dataset = file.createDataSet(array_name, datatype, dataspace);

  // Write to dataset using default memory space
  dataset.write(array.data(), datatype, dataspace);
}
