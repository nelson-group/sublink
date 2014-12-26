#pragma once

/** @file ReadSubfindHDF5.hpp
 * @brief Utilities to read Subfind HDF5 snapshot files.
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

#include <vector>
#include <string>
#include <sstream>
#include <cassert>

#include "H5Cpp_wrapper.hpp"

/** @namespace subfind
 * @brief Namespace containing functions to read Subfind output files.
 */
namespace subfind {

// Number of particle types in the simulation
static constexpr int num_parttypes = 6;

/** @brief Function to read a scalar attribute from a Subfind output file.
 *
 * @tparam T Type of the attribute.
 * @param[in] basedir Directory containing the Subfind output files.
 * @param[in] snapnum Snapshot number.
 * @param[in] attr_name Name of the attribute.
 * @param[in] filenum File number.
 * @return The attribute value.
 */
template <class T>
T get_scalar_attribute(const std::string& basedir, const int16_t snapnum,
    const std::string& attr_name, const int32_t filenum = 0) {

  // Create filename
  std::stringstream tmp_stream;
  tmp_stream << basedir << "/groups_" <<
      std::setfill('0') << std::setw(3) << snapnum << "/fof_subhalo_tab_" <<
      std::setfill('0') << std::setw(3) << snapnum;
  tmp_stream << "." << filenum << ".hdf5";
  std::string file_name = tmp_stream.str();

  // Open attribute from file header
  auto file = new H5::H5File(file_name, H5F_ACC_RDONLY );
  auto header_group = new H5::Group(file->openGroup("Header"));
  auto myattr = new H5::Attribute(header_group->openAttribute(attr_name));

  // Check that type sizes match (necessary but not sufficient).
  auto dt = myattr->getDataType();
  assert(dt.getSize() == sizeof(T));

  // Read attribute
  T retval;
  myattr->read(myattr->getDataType(), &retval);

  // Close file and release resources
  file->close();
  delete myattr;
  delete header_group;
  delete file;

  return retval;
}

/** @brief Function to read a one-dimensional dataset (a.k.a. block)
 *         or a dataset "column" from a single file.
 *
 * @tparam T Type of the elements in the dataset.
 * @param[in] file_name Path to the input file.
 * @param[in] group_name The HDF5 group name, i.e., @a Group or @a Subhalo.
 * @param[in] block_name Name of the dataset.
 * @param[in] parttype The particle type; -1 for one-dimensional datasets.
 *
 * @pre Dataset is one- or two-dimensional.
 * @pre If the dataset is two-dimensional, the number of columns
 *      corresponds to the number of particle types.
 *
 * @return A vector with the dataset values.
 */
template <class T>
std::vector<T> read_block_single_file(const std::string& file_name,
    const std::string& group_name, const std::string& block_name,
    const int parttype = -1) {

  // Open the specified dataset.
  auto file = new H5::H5File(file_name, H5F_ACC_RDONLY );
  auto group = new H5::Group(file->openGroup(group_name));
  auto dataset = new H5::DataSet(group->openDataSet(block_name));

  // Get dimensions of the dataset
  H5::DataSpace file_space = dataset->getSpace();
  const unsigned int rank = file_space.getSimpleExtentNdims();
  hsize_t dims_out[rank];
  file_space.getSimpleExtentDims(dims_out, NULL);

  // Check that dataspace makes sense.
  if (parttype == -1)
    assert(rank == 1);
  else {
    assert((rank == 2) && (dims_out[1] == num_parttypes));
    // Define hyperslab for column of interest.
    hsize_t offset[2] = {0, parttype};     // hyperslab offset in file
    hsize_t count[2]  = {dims_out[0], 1};  // hyperslab size in file
    file_space.selectHyperslab(H5S_SELECT_SET, count, offset);
  }

  // Create output dataspace (in memory)
  hsize_t dimsm[1] = {dims_out[0]};
  H5::DataSpace mem_space(1, dimsm);

  // Read data to vector.
  std::vector<T> retval(dimsm[0]);
  dataset->read(retval.data(), dataset->getDataType(), mem_space, file_space);

  // Close file and release resources
  file->close();
  delete dataset;
  delete group;
  delete file;

  return retval;
}

/** @brief Function to read a one-dimensional dataset (a.k.a. block)
 *         or a dataset "column."
 *
 * @tparam T Type of the elements in the dataset.
 * @param[in] basedir Directory containing the Subfind output files.
 * @param[in] snapnum Snapshot number.
 * @param[in] group_name The HDF5 group name, i.e., @a Group or @a Subhalo.
 * @param[in] block_name Name of the dataset.
 * @param[in] parttype The particle type; -1 for one-dimensional datasets.
 *
 * @pre (group_name == "Group") || (group_name == "Subhalo")
 *
 * @return A vector with the dataset values.
 */
template <class T>
std::vector<T> read_block(const std::string& basedir,
    const int16_t snapnum, const std::string& group_name,
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
    if (len_thisfile <= 0) {
      std::cout << "Group file " << filenum << " from snapshot " << snapnum <<
          " is empty.\n";
      continue;
    }

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
