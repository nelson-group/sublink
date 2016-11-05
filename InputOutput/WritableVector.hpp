#pragma once

/** @file WritableVector.hpp
 * @brief Defines the WritableVector class, which is convenient for writing
 *        data to HDF5 files.
 *
 * @author Vicente Rodriguez-Gomez (vrg@jhu.edu)
 */

#include <vector>
#include <string>
#include <cassert>

#include "H5Cpp_wrapper.hpp"

/** @class WritableVector
 * @brief Vector wrapper that makes the process of writing to HDF5 files
 *        more straightforward.
 *
 * @note To be implemented later. Not currently used.
 */
template <typename T>
class WritableVector {

public:
  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  /** Default constructor. Creates invalid WritableVector. */
  WritableVector() : v_(), name_(), dt_() {
  }

  /** Constructor. */
  WritableVector(const std::vector<T>& v, const std::string& name,
      const H5::DataType& dt) : v_(v), name_(name), dt_(dt) {
  }

  std::vector<T>& vector() {
    return v_;
  }

  const std::vector<T>& vector() const {
    return v_;
  }

  const std::string& name() const {
    return name_;
  }

  const H5::DataType& datatype() const {
    return dt_;
  }

private:
  std::vector<T> v_;
  std::string name_;
  H5::DataType dt_;
};
