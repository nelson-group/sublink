#pragma once
/** @file ReadArepoHDF5.hpp
 * @brief Utilities to read Arepo HDF5 snapshot files.
 *
 * Usage example 1: read particle IDs and masses of star particles
 * @code{.cpp}
 * std::vector<uint64_t> all_ids;
 * std::vector<float> all_masses;
 * std::string snapname = "snapdir_135/snap_135";
 * const int parttype = 4;
 * all_ids = arepo::read_block<uint64_t>(snapname, "ParticleIDs", parttype);
 * all_masses = arepo::read_block<float>(snapname, "Masses", parttype);
 * @endcode
 *
 * Usage example 2: read positions of DM particles (define datastruct)
 * @code{.cpp}
 * #include <array>
 * using PositionType = std::array<float, 3>;
 * void main() {
 *   std::string snapname = "snapdir_135/snap_135";
 *   const int parttype = 1;
 *   std::vector<PositionType> pos;
 *   pos = arepo::read_block<PositionType>(snapname, "Coordinates", parttype);
 * }
 * @endcode
 *
 * @author Vicente Rodriguez-Gomez (vrg@jhu.edu)
 */

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>  // setfill, setw, setprecision
#include <cassert>

#include "H5Cpp_wrapper.hpp"
#include "../Util/TreeTypes.hpp"

/** @namespace arepo
 * @brief Namespace containing functions to read Arepo snapshot files.
 */
namespace arepo {

/** @brief Function to read a scalar attribute from the snapshot header.
 *
 * @tparam T Type of the attribute.
 * @param[in] file_name Path to the input file.
 * @param[in] attr_name Name of the attribute.
 * @return The attribute value.
 */
template <typename T>
T get_scalar_attribute(const std::string& file_name,
    const std::string& attr_name) {

  // Open attribute from snapshot header
  auto file = H5::H5File(file_name, H5F_ACC_RDONLY );
  auto header_group = H5::Group(file.openGroup("Header"));
  auto myattr = H5::Attribute(header_group.openAttribute(attr_name));

  // Check that type sizes match (necessary but not sufficient).
  auto dt = myattr.getDataType();
  assert(dt.getSize() == sizeof(T));

  // Read attribute
  T retval;
  myattr.read(myattr.getDataType(), &retval);

  file.close();

  return retval;
}

/** @brief Function to read an attribute array from the snapshot header.
 *
 * @tparam T Type of the elements of the attribute array.
 * @param[in] file_name Path to the input file.
 * @param[in] attr_name Name of the attribute array.
 * @return A vector with the attribute array values.
 */
template <typename T>
std::vector<T> get_vector_attribute(const std::string& file_name,
    const std::string& attr_name) {
  // Open attribute from snapshot header
  auto file = H5::H5File(file_name, H5F_ACC_RDONLY );
  auto header_group = H5::Group(file.openGroup("Header"));
  auto myattr = H5::Attribute(header_group.openAttribute(attr_name));

  // Check that type sizes match (necessary but not sufficient).
  auto dt = myattr.getDataType();
  assert(dt.getSize() == sizeof(T));

  // Deduce length of attribute.
  auto space = myattr.getSpace();
  assert(space.getSimpleExtentNdims() == 1);
  auto attr_len = space.getSimpleExtentNpoints();

  // Read attribute values directly into vector.
  std::vector<T> retval(attr_len);
  myattr.read(myattr.getDataType(), retval.data());

  file.close();

  return retval;
}

/** @brief Get the datatype size of a given dataset (a.k.a. block).
 * @param[in] basedir Directory containing the snapshot files.
 * @param[in] snapnum Snapshot number.
 * @param[in] block_name Name of the dataset.
 * @param[in] parttype The particle type.
 *
 * @return The datatype size of the dataset.
 */
std::size_t get_datatype_size(const std::string& basedir,
    const int16_t snapnum, const std::string& block_name,
    const int parttype) {

  // Group name corresponding to particle type.
  std::stringstream ss;
  ss << "/PartType" << parttype;
  std::string parttype_str = ss.str();

  // Snapshot filename without the file number
  ss.str("");
  ss << basedir << "/snapdir_" <<
      std::setfill('0') << std::setw(3) << snapnum << "/snap_" <<
      std::setfill('0') << std::setw(3) << snapnum;
  std::string file_name_base = ss.str();

  // Name of first file (just to read header)
  ss.str("");
  ss << file_name_base << ".0.hdf5";
  std::string file_name = ss.str();

  // Iterate over files until we find one that contains the specified dataset.
  auto nfiles = get_scalar_attribute<int32_t>(file_name, "NumFilesPerSnapshot");
  for (int32_t filenum=0; filenum < nfiles; filenum++) {
    // Open current file
    ss.str("");
    ss << file_name_base << "." << filenum << ".hdf5";
    file_name = ss.str();
    auto file = H5::H5File(file_name, H5F_ACC_RDONLY);
    // Only proceed if group exists.
    if (!H5Lexists(file.getId(), parttype_str.data(), H5P_DEFAULT)) {
      file.close();
      continue;
    }
    // Return datatype size.
    auto group = H5::Group(file.openGroup(parttype_str));
    auto dataset = H5::DataSet(group.openDataSet("Coordinates"));
    auto dt = dataset.getDataType();
    file.close();
    return dt.getSize();
  }
  std::cerr << "ERROR: " << parttype_str << "/" << block_name << " not found!\n";
  assert(false);
  return 0;
}


/** @brief Function to read a dataset (a.k.a. block)
 * from a single file.
 *
 * @tparam T Type of the elements in the dataset.
 * @param[in] file_name Path to the input file.
 * @param[in] block_name Name of the dataset.
 * @param[in] parttype The particle type.
 * @param[in] read_num If nonzero, read -at most- this many elements.
 * @return A vector with the dataset values.
 *
 * @note The return value is a vector (with an appropriately chosen type).
 *       When reading from a dataset with ndims > 1, the size() of the
 *       output vector equals the size of the first dimension of the dataset.
 */
template <typename T>
std::vector<T> read_block_single_file(const std::string& file_name,
    const std::string& block_name, const int parttype, const uint64_t read_num = 0) {

  // Create group name corresponding to particle type.
  std::stringstream ss;
  ss << "/PartType" << parttype;
  std::string parttype_str = ss.str();

  // Only proceed if group exists (have a peek at the file).
  auto tmp_file = H5::H5File(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (H5Lexists(tmp_file.getId(), parttype_str.data(), H5P_DEFAULT) == false) {
    tmp_file.close();
    return std::vector<T>();
  }
  tmp_file.close();

  // Open the specified dataset.
  auto file = H5::H5File(file_name, H5F_ACC_RDONLY );
  auto group = H5::Group(file.openGroup(parttype_str));
  auto dataset = H5::DataSet(group.openDataSet(block_name));

  // Get dimensions of the dataset
  H5::DataSpace file_space = dataset.getSpace();
  const unsigned int file_rank = file_space.getSimpleExtentNdims();
  hsize_t file_dims[file_rank];
  file_space.getSimpleExtentDims(file_dims, NULL);

  // If read_num != 0, then partial read only, of -at most- this number of (1D) elements
  if ((read_num > 0) && (read_num < file_dims[0])) {
    assert(file_rank == 1);
    file_dims[0] = read_num;
    hsize_t count[1] = {read_num};
    hsize_t offset[1] = {0};
    file_space.selectHyperslab( H5S_SELECT_SET, count, offset );
  }

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

/** @brief Function to read a dataset (a.k.a. block) from the snapshot files.
 *
 * @tparam T Type of the elements in the datasets.
 * @param[in] basedir Directory containing the snapshot files.
 * @param[in] snapnum Snapshot number.
 * @param[in] block_name Name of the dataset.
 * @param[in] parttype The particle type.
 * @param[in] read_num If nonzero, only read this many total elements (1D).
 * @return A vector with the dataset values.
 *
 * @note The return value is a vector (with an appropriately chosen type).
 *       When reading from a dataset with ndims > 1, the size() of the
 *       output vector equals the size of the first dimension of the
 *       (concatenated) dataset.
 */
template <typename T>
std::vector<T> read_block(const std::string& basedir,
    const int16_t snapnum, const std::string& block_name,
    const int parttype, const uint64_t read_num = 0) {

  // Snapshot filename without the file number
  std::stringstream tmp_stream;
  tmp_stream << basedir << "/snapdir_" <<
      std::setfill('0') << std::setw(3) << snapnum << "/snap_" <<
      std::setfill('0') << std::setw(3) << snapnum;
  std::string file_name_base = tmp_stream.str();

  // Name of first file (just to read header)
  std::stringstream sstream;
  sstream << file_name_base << ".0.hdf5";
  std::string file_name = sstream.str();

  // Get total number of particles
  auto npart_total_vect_lowword = get_vector_attribute<uint32_t>(file_name,
      "NumPart_Total");
  auto npart_total_vect_highword = get_vector_attribute<uint32_t>(file_name,
      "NumPart_Total_HighWord");
  part_id_type npart_total = (int64_t)npart_total_vect_lowword[parttype] | 
                             ((int64_t)npart_total_vect_highword[parttype] << 32);

  // Modify to selective global read
  if(read_num > 0)
    npart_total = read_num;

  // Allocate memory for output vector
  std::vector<T> data_total;
  data_total.reserve(npart_total);

  // Load data from each file and append to big vector
  auto nfiles = get_scalar_attribute<int32_t>(file_name, "NumFilesPerSnapshot");
  for (int32_t filenum=0; filenum < nfiles; filenum++) {
    // Load data from current file.
    sstream.str("");
    sstream << file_name_base << "." << filenum << ".hdf5";
    file_name = sstream.str();
    auto data_thisfile = read_block_single_file<T>(file_name, block_name, parttype, npart_total);

    // Get number of particles according to file header (just to check)
    auto npart_thisfile_vect = get_vector_attribute<int32_t>(file_name,
        "NumPart_ThisFile");
    part_id_type npart_thisfile = npart_thisfile_vect[parttype];
    if ((data_thisfile.size() != npart_thisfile) && (read_num == 0))
      std::cerr << "BAD: number of particles in file " << filenum <<
        " does not match with header.\n";

    // Concatenate vectors.
    data_total.insert(data_total.end(), data_thisfile.begin(), data_thisfile.end());

    npart_total -= data_thisfile.size(); // number left to read
    if(npart_total == 0)
      break;
  }

  // Check total number of particles with header.
  if ((data_total.size() != npart_total) && (npart_total > 0))
    std::cerr << "BAD: total number of particles does not match with header ["
              << data_total.size() << " vs " << npart_total << "].\n";

  return data_total;
}

}  // end namespace arepo
