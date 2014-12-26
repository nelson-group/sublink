#pragma once

/** @file ReadWriteHDF5.hpp
 * @brief General-purpose utilities to read from and write to HDF5 files.
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

#include <vector>
#include <string>
#include <cassert>

#include "H5Cpp_wrapper.hpp"

/** @brief Function to read a one-dimensional dataset (a.k.a. block)
 * from a single HDF5 file.
 *
 * @tparam T Type of the elements in the dataset.
 * @param[in] file_name Path to the input file.
 * @param[in] block_name Name of the dataset.
 * @return A vector with the dataset values.
 */
template <typename T>
std::vector<T> read_dataset(const std::string& file_name,
    const std::string& block_name) {

  // Open dataset
  H5::H5File file(file_name, H5F_ACC_RDONLY );
  H5::DataSet dataset(file.openDataSet(block_name));

  // Get dimensions of the dataset
  H5::DataSpace file_space = dataset.getSpace();
  const unsigned int rank = file_space.getSimpleExtentNdims();
  hsize_t dims_out[rank];
  file_space.getSimpleExtentDims(dims_out, NULL);

  // Dataspace in memory is the same as in the HDF5 file.
  hsize_t dimsm[rank];
  for (hsize_t k = 0; k < rank; k++) dimsm[k] = dims_out[k];
  H5::DataSpace mem_space(rank, dimsm);

  // Read data
  std::vector<T> retval(dimsm[0]);  // Output vector is 1D.
//  dataset.read(&retval[0], dataset.getDataType(), mem_space, file_space);
  dataset.read(retval.data(), dataset.getDataType(), mem_space, file_space);

  file.close();
  return retval;
}

/** @brief Function to open an HDF5 file in "append" mode and add
 *         a new one-dimensional array.
 *
 *
 * @note Byte order is little-endian by default.
 */
template <typename T>
void add_array(const std::string& writefilename,
    const std::vector<T>& array,
    const std::string& array_name,
    H5::DataType datatype) {

    // Open file with read-write access (create if it does not exist).
    H5::H5File* file = new H5::H5File(writefilename, H5F_ACC_TRUNC);

    // Only proceed if array is non-empty
    if (array.size() == 0) {
      file->close();
      delete file;
      return;
    }

    // Define (one-dimensional) dataspace
    hsize_t dimsf[1];  // dataset dimensions
    dimsf[0] = array.size();
    H5::DataSpace dataspace(1, dimsf);  // rank 1

//    // Define datatypes (little-endian) for int8, int16, etc.
//    H5::IntType datatype_int8(H5::PredType::NATIVE_INT8);
//    H5::IntType datatype_int16(H5::PredType::NATIVE_INT16);
//    H5::IntType datatype_int32(H5::PredType::NATIVE_INT32);
//    H5::IntType datatype_uint64(H5::PredType::NATIVE_UINT64);
//    datatype_int8.setOrder(H5T_ORDER_LE);
//    datatype_int16.setOrder(H5T_ORDER_LE);
//    datatype_int32.setOrder(H5T_ORDER_LE);
//    datatype_uint64.setOrder(H5T_ORDER_LE);

    // Create dataset
    H5::DataSet dataset = file->createDataSet(array_name, datatype, dataspace);

    // Write to dataset using default memory space
    dataset.write(array.data(), datatype, dataspace);

    // Clear
    file->close();
    delete file;
  }



}



