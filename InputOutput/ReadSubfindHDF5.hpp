#pragma once
/** @file ReadSubfindHDF5.hpp
 * @brief Utilities to read Subfind HDF5 snapshot files.
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>  // setfill, setw, setprecision
#include <cassert>

#include "H5Cpp_wrapper.hpp"
#include "../Util/TreeTypes.hpp"

/** @namespace subfind
 * @brief Namespace containing functions to read Subfind output files.
 */
namespace subfind {

/** @brief Function to read a scalar attribute from a Subfind output file.
 *
 * @tparam T Type of the attribute.
 * @param[in] basedir Directory containing the Subfind output files.
 * @param[in] snapnum Snapshot number.
 * @param[in] attr_name Name of the attribute.
 * @param[in] filenum File number.
 * @return The attribute value.
 */
template <typename T>
T get_scalar_attribute(const std::string& basedir, const snapnum_type snapnum,
    const std::string& attr_name, const int32_t filenum = 0) {

  // Create filename
  std::stringstream tmp_stream;
  tmp_stream << basedir << "/groups_" <<
      std::setfill('0') << std::setw(3) << snapnum << "/fof_subhalo_tab_" <<
      std::setfill('0') << std::setw(3) << snapnum;
  tmp_stream << "." << filenum << ".hdf5";
  std::string file_name = tmp_stream.str();

  // Open attribute from file header
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

/** @brief Function to read a dataset (a.k.a. block)
 *         or a dataset "column" from a single file.
 *
 * @tparam T Type of the elements in the dataset.
 * @param[in] file_name Path to the input file.
 * @param[in] group_name The HDF5 group name, i.e., @a Group or @a Subhalo.
 * @param[in] block_name Name of the dataset.
 * @param[in] parttype The 0-indexed column number, if any (-1 otherwise).
 *
 * @pre Dataset is one- or two-dimensional.
 * @return A vector with the dataset values.
 */
template <typename T>
std::vector<T> read_block_single_file(const std::string& file_name,
    const std::string& group_name, const std::string& block_name,
    const int parttype = -1) {

  // Open the specified dataset.
  auto file = H5::H5File(file_name, H5F_ACC_RDONLY );
  auto group = H5::Group(file.openGroup(group_name));
  auto dataset = H5::DataSet(group.openDataSet(block_name));

  // Get dimensions of the dataset
  H5::DataSpace file_space = dataset.getSpace();
  const unsigned int rank = file_space.getSimpleExtentNdims();
  hsize_t dims_out[rank];
  file_space.getSimpleExtentDims(dims_out, NULL);

  // When parttype == -1 (default), dataspace in memory is the same
  // as in the HDF5 file.
  hsize_t dimsm[rank];
  for (hsize_t k = 0; k < rank; ++k) dimsm[k] = dims_out[k];
  H5::DataSpace mem_space(rank, dimsm);

  // Output vector is 1D.
  std::vector<T> retval(dimsm[0]);

  // If parttype != -1, we are only interested in a single column.
  if (parttype != -1) {
    assert((rank == 2) && (parttype >= 0));
    // Define hyperslab (offset and size) for column of interest.
    hsize_t offset[2] = {0, static_cast<hsize_t>(parttype)};
    hsize_t count[2]  = {dims_out[0], 1};
    file_space.selectHyperslab(H5S_SELECT_SET, count, offset);
    // The dataspace becomes 1D
    hsize_t dimsm_1d[1] = {dims_out[0]};
    mem_space.setExtentSimple(1, dimsm_1d);
  }

  // Read data
  dataset.read(retval.data(), dataset.getDataType(), mem_space, file_space);

  // Close file and release resources
  file.close();

  return retval;
}

/** @brief Function to read a dataset (a.k.a. block)
 *         or a dataset "column."
 *
 * @tparam T Type of the elements in the dataset.
 * @param[in] basedir Directory containing the Subfind output files.
 * @param[in] snapnum Snapshot number.
 * @param[in] group_name The HDF5 group name, i.e., @a Group or @a Subhalo.
 * @param[in] block_name Name of the dataset.
 * @param[in] parttype The 0-indexed column number, if any (-1 otherwise).
 *
 * @pre (group_name == "Group") || (group_name == "Subhalo")
 *
 * @return A vector with the dataset values.
 */
template <typename T>
std::vector<T> read_block(const std::string& basedir,
    const snapnum_type snapnum, const std::string& group_name,
    const std::string& block_name, const int parttype = -1) {

  // Check precondition
  assert((group_name == "Group") || (group_name == "Subhalo"));

  // Create initial filename
  std::stringstream tmp_stream;
  tmp_stream << basedir << "/groups_" <<
      std::setfill('0') << std::setw(3) << snapnum << "/fof_subhalo_tab_" <<
      std::setfill('0') << std::setw(3) << snapnum;
  std::string file_name_base = tmp_stream.str();
  tmp_stream << ".0.hdf5";
  std::string file_name = tmp_stream.str();

  // Get total number of objects
  std::string attr_name;
  if (group_name == "Group")
    attr_name = "Ngroups_Total";
  else
    attr_name = "Nsubgroups_Total";
  auto len_total = get_scalar_attribute<int32_t>(basedir, snapnum, attr_name);

  // Allocate memory for result
  std::vector<T> data_total;
  data_total.reserve(len_total);

  // Iterate over files to read data
  auto nfiles = get_scalar_attribute<int32_t>(basedir, snapnum, "NumFiles");
  for (int32_t filenum = 0; filenum < nfiles; ++filenum) {

    // Get number of objects in this file
    if (group_name == "Group")
      attr_name = "Ngroups_ThisFile";
    else
      attr_name = "Nsubgroups_ThisFile";
    auto len_thisfile = get_scalar_attribute<int32_t>(basedir, snapnum,
        attr_name, filenum);

    // Skip empty group files
    if (len_thisfile <= 0)
      continue;

    // Create filename and read data
    tmp_stream.str("");
    tmp_stream << file_name_base << "." << filenum << ".hdf5";
    file_name = tmp_stream.str();
    std::vector<T> data_thisfile = read_block_single_file<T>(file_name,
        group_name, block_name, parttype);

    // Sanity check
    if (data_thisfile.size() != static_cast<uint32_t>(len_thisfile))
      std::cerr << "BAD: number of objects in file " << filenum <<
                   " does not match with header.\n";

    // Append data to output array
    data_total.insert(data_total.end(), data_thisfile.begin(),
        data_thisfile.end());
  }

  // Check total number of particles
  if (data_total.size() != static_cast<uint32_t>(len_total))
    std::cerr << "BAD: total number of objects does not match with header.\n";

  return data_total;
}

}  // end namespace subfind
