#pragma once

/** @file ReadArepoHDF5.hpp
 * @brief Utilities to read Arepo HDF5 snapshot files.
 *
 * Usage example 1: read particle IDs and masses of star particles
 * std::vector<uint64_t> all_ids;
 * std::vector<float> all_masses;
 * std::string snapname = "snapdir_135/snap_135";
 * const int parttype = 4;
 * all_ids = arepo_read_block<uint64_t>(snapname, "ParticleIDs", parttype);
 * all_masses = arepo_read_block<float>(snapname, "Masses", parttype);
 *
 * Usage example 2: read positions of DM particles (define datastruct)
 * struct pos_struct{
 *    float x;
 *    float y;
 *    float z;
 *  };
 *  void main() {
 *    std::string snapname = "snapdir_135/snap_135";
 *    const int parttype = 1;
 *    std::vector<pos_struct> pos;
 *    pos = arepo_read_block<pos_struct>(snapname, "Coordinates", parttype);
 *  }
 *
 * Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

#include <vector>
#include <string>
#include <sstream>


// Note that "inttypes.h" is included by H5Cpp.h,
// so hopefully there is no conflict.
#include "H5Cpp_wrapper.hpp"

/** Function to read a scalar attribute from the header
 * of a snapshot file. */
template <class T>
T arepo_get_attribute(std::string file_name, std::string attr_name) {
  // Open attribute from header
  H5::H5File* file = new H5::H5File(file_name, H5F_ACC_RDONLY );
  H5::Group* header_group = new H5::Group(file->openGroup("Header"));
  H5::Attribute* myattr = new H5::Attribute (header_group->openAttribute(attr_name));

  // Read attribute
  T attr_value;
  myattr->read(myattr->getDataType(), &attr_value);

  // Close file and release resources
  file->close();
  delete myattr;
  delete header_group;
  delete file;

  return attr_value;
}

/** Function to read an attribute with length 6 from the header
 * of a snapshot file.
 *
 * @note: Only works with integer types.
 */
template <class T>
std::vector<T> arepo_get_attribute_array(std::string file_name, std::string attr_name) {
  // Open attribute from header
  H5::H5File* file = new H5::H5File(file_name, H5F_ACC_RDONLY );
  H5::Group* header_group = new H5::Group(file->openGroup("Header"));
  H5::Attribute* myattr = new H5::Attribute (header_group->openAttribute(attr_name));

  // Read attribute into dynamic array, and then convert into vector
  // (in this case, H5 refuses to read directly into the vector)
  T* tmp_array = new T;
  myattr->read(myattr->getDataType(), tmp_array);

  int num_parttypes = 6; // there are always 6 types of particles
  std::vector<T> npart_total( tmp_array, tmp_array + num_parttypes);

  // Close file and release resources
  file->close();
  delete tmp_array;
  delete myattr;
  delete header_group;
  delete file;

  return npart_total;
}

/** Function to read a one-dimensional block from a single file. */
template <class T>
std::vector<T> arepo_read_block_single_file(const std::string file_name,
    const std::string block_name, const int parttype) {

  // Define particle type
  std::stringstream sstream;
  sstream << "PartType" << parttype;

  // ------------------------------------------------------
  // Fix (2013/09/06): only proceed if group exists.
  // Need the group path (e.g. "/PartType0") as a char array:
  std::string str = "/" + sstream.str();
  const char *cstr = str.c_str();
  // Have a peek at the file to see if group exists
  H5::H5File tmp_file;
  tmp_file = H5::H5File(file_name, H5F_ACC_RDONLY, H5P_DEFAULT );
  if (H5Lexists(tmp_file.getId(), cstr, H5P_DEFAULT) == false) {
    std::cout << "Group " << sstream.str() << " does not exist. Skipping...\n";
    tmp_file.close();
    // Return empty vector
    return std::vector<T>();
  }
  tmp_file.close();
  // ------------------------------------------------------

  // Open data set with given block name
  H5::H5File* file = new H5::H5File(file_name, H5F_ACC_RDONLY );
  H5::Group* group = new H5::Group(file->openGroup(sstream.str()));
  H5::DataSet* dataset = new H5::DataSet(group->openDataSet(block_name));

  // Get dataspace of the data set
  H5::DataSpace dataspace = dataset->getSpace();

  const unsigned int rank = dataspace.getSimpleExtentNdims();
  hsize_t *dims_out = new hsize_t[rank];
  dataspace.getSimpleExtentDims(dims_out, NULL);

  // Read data from file into data_out (in memory)
  hsize_t *dimsm = new hsize_t[rank];
  for (unsigned int k = 0; k < rank; k++) {
    dimsm[k] = dims_out[k];
  }

  H5::DataSpace memspace(rank, dimsm);
  std::vector<T> data_out(dimsm[0]);
  dataset->read( &(data_out[0]), dataset->getDataType(), memspace, dataspace);

  // Close file and release resources
  file->close();
  delete[] dims_out;
  delete[] dimsm;
  delete dataset;
  delete group;
  delete file;

  return data_out;
}

/** Function to read a one-dimensional block. */
template <class T>
std::vector<T> arepo_read_block(const std::string file_name_base,
    const std::string block_name, const int parttype) {
  // NOTE: In most simulations, particle IDs have type uint64.
  // However, the attributes NumPart_Total and NumPart_ThisFile have
  // types uint32 and int32, respectively, even for the larger simulations...

  // Initial file name (to read header)
  std::stringstream sstream;
  sstream << file_name_base << ".0.hdf5";
  std::string file_name = sstream.str();

  // Get total number of particles
  std::vector<unsigned int> npart_total_vect = arepo_get_attribute_array<unsigned int>(
      file_name, "NumPart_Total");
  unsigned int npart_total = npart_total_vect[parttype];

  // Allocate memory for result (huge vector)
  std::vector<T> data_total;
  data_total.reserve(npart_total);

  // Get number of files
  int nfiles = arepo_get_attribute<int>(file_name, "NumFilesPerSnapshot");

  // Load data from each file and append to big vector
  std::vector<T> data_thisfile;
  std::vector<int> npart_thisfile_vect; // just to check
  for (int filenum=0; filenum < nfiles; filenum++) {
    // create filename
    sstream.str("");
    sstream << file_name_base << "." << filenum << ".hdf5";
    file_name = sstream.str();
    // load IDs from this file
    data_thisfile = arepo_read_block_single_file<T>(file_name, block_name, parttype);
    // load number of particles from file header (just to check)
    npart_thisfile_vect = arepo_get_attribute_array<int>(file_name, "NumPart_ThisFile");
    if (data_thisfile.size() != (unsigned int)npart_thisfile_vect[parttype]) {
      std::cerr << "BAD: number of particles in file " << filenum << " do not match.\n";
    }
    // Concatenate vectors
    data_total.insert(data_total.end(), data_thisfile.begin(), data_thisfile.end());
  }
  // Check total number of particles
  if (data_total.size() != npart_total) {
    std::cerr << "BAD: total number of particles does not match with header.\n";
  }

  return data_total;
}


