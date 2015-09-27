#pragma once
/** @file AllGalaxies.hpp
 * @brief Define the AllGalaxies and Galaxy classes.
 *
 * @author Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)
 */

#include <tuple>  // std::tie
#include <cassert>
#include "Point.hpp"

/** Type of indexes and sizes. */
typedef uint32_t size_type;

/** @class AllGalaxies
 * @brief Define container for all galaxies in the simulation snapshot.
 *
 * AF(*this) = <G>, where:
 *   G = [g_0, g_1, ..., g_{n-1}], n = num_galaxies().
 */
template <typename V = char>
class AllGalaxies {
private:
  // Predeclare internal type for galaxies.
  struct internal_galaxy;

public:
  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  /** Type of this container. */
  typedef AllGalaxies all_galaxies_type;

  /** Predeclaration of Galaxy type. */
  class Galaxy;
  /** Synonym for Galaxy (following STL conventions). */
  typedef Galaxy galaxy_type;

  /** Type of galaxy iterators. */
  class GalaxyIterator;
  /** Synonym for GalaxyIterator */
  typedef GalaxyIterator galaxy_iterator;

  /** Galaxy value type. */
  typedef V galaxy_value_type;

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  /** Construct an empty container. */
  AllGalaxies()
      : galaxies_() {
  }

  /** Default destructor */
  ~AllGalaxies() = default;

  /////////////
  // General //
  /////////////

  /** Return the total number of galaxies.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return galaxies_.size();
  }

  /** Remove all galaxies.
   * @post num_galaxies() == 0
   *
   * Invalidates all outstanding Galaxy objects.
   */
  void clear() {
    galaxies_.clear();
  }

  ///////////////
  // GALAXIES //
  ///////////////

  /** @class AllGalaxies::Galaxy
   * @brief Class to represent galaxies.
   *
   * AF(p) = g_i = <g_i.position_, g_i.value>, where
   *   g_i is the ith galaxy in AF(*p.ag_).
   */
  class Galaxy {
  public:

    /** Construct invalid Galaxy. */
    Galaxy() : ag_(nullptr), idx_(-1) {
    }

    /** Return this Galaxy's position. */
    Point& position() {
      return fetch().point_;
    }
    /** Return this Galaxy's position, read-only. */
    const Point& position() const {
      return fetch().point_;
    }
    /** Return this Galaxy's index, a number in the range
     * [0, ag_->num_galaxies()).
     */
    size_type index() const {
      assert(is_valid());
      return idx_;
    }
    /** Return value associated with this Galaxy. */
    galaxy_value_type& value() {
      return fetch().value_;
    }
    /** Return value associated with this Galaxy, read-only. */
    const galaxy_value_type& value() const {
      return fetch().value_;
    }

    /** Test whether this Galaxy and @a x are equal. */
    bool operator==(const Galaxy& x) const {
      return ((ag_ == x.ag_) && (idx_ == x.idx_));
    }
    /** Test whether this Galaxy is less than @a x in the global order.
     *
     * The ordering relation must obey trichotomy:
     * For any two Galaxies x and y, exactly one of x == y,
     * x < y, and y < x is true.
     */
    bool operator<(const Galaxy& x) const {
      // Compare ag_, then idx_ (pro tip accepted!)
      return std::tie(ag_, idx_) < std::tie(x.ag_, x.idx_);
    }

  private:
    friend class AllGalaxies;
    // Pointer back to the AllGalaxies container.
    all_galaxies_type* ag_;
    // Index of this galaxy.
    size_type idx_;

    /** Private constructor. */
    Galaxy(const all_galaxies_type* ag, size_type idx)
        : ag_(const_cast<all_galaxies_type*>(ag)), idx_(idx) {
    }
    /** Helper method to determine whether this Galaxy is valid. */
    bool is_valid() const {
      return (ag_ != nullptr) && (idx_ < ag_->size());
    }
    /** Helper method to return a reference to the internal galaxy.
     * @pre This galaxy is valid.
     */
    internal_galaxy& fetch() const {
      assert(is_valid());
      return ag_->galaxies_[idx_];
    }
  };

  /** Synonym for size(). */
  size_type num_galaxies() const {
    return size();
  }

  /** Add a Galaxy to this container, returning the added Galaxy.
   * @param[in] position The new Galaxy's position
   * @param[in] value The new Galaxy's value (optional)
   * @post new size() == old size() + 1
   * @post result.index() == old size()
   *
   * Complexity: O(1) amortized operations.
   */
  Galaxy add_galaxy(const Point& position,
      const galaxy_value_type& value = galaxy_value_type()) {
    galaxies_.emplace_back(internal_galaxy(position, value));
    return Galaxy(this, galaxies_.size()-1);
  }

  ///////////////
  // Iterators //
  ///////////////

  /** @class AllGalaxies::GalaxyIterator
   * @brief Iterator class for galaxies. A forward iterator. */
  class GalaxyIterator : private totally_ordered<GalaxyIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Galaxy value_type;
    /** Type of pointers to elements. */
    typedef Galaxy* pointer;
    /** Type of references to elements. */
    typedef Galaxy& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid GalaxyIterator. */
    GalaxyIterator() : ag_(nullptr), idx_(-1) {
    }
    /** Method to dereference a galaxy iterator.
     * @pre This is a valid GalaxyIterator.
     */
    Galaxy operator*() const {
      assert(is_valid());
      return Galaxy(ag_, idx_);
    }
    /** Method to increment a galaxy iterator.
     * @pre This is a valid GalaxyIterator.
     */
    GalaxyIterator& operator++() {
      assert(is_valid());
      ++idx_;
      return *this;
    }
    /** Method to compare two galaxy iterators. */
    bool operator==(const GalaxyIterator& x) const {
      return ((idx_ == x.idx_) && ag_ == x.ag_);
    }

   private:
    friend class AllGalaxies;
    // Pointer back to the AllGalaxies container.
    all_galaxies_type* ag_;
    // Index of this galaxy.
    size_type idx_;

    /** Private constructor. */
    GalaxyIterator(const all_galaxies_type* ag, size_type idx)
        : ag_(const_cast<all_galaxies_type*>(ag)), idx_(idx) {
    }
    /** Helper method to determine whether this GalaxyIterator is valid.
     */
    bool is_valid() const {
      return (ag_ != nullptr) && (idx_ < ag_->size());
    }
  };

  /** Return an iterator pointing to the first galaxy. */
  galaxy_iterator begin() const {
    assert(num_galaxies() > 0);
    return GalaxyIterator(this, 0);
  }
  /** Return an iterator pointing to one-past-the-last galaxy. */
  galaxy_iterator end() const {
    return GalaxyIterator(this, num_galaxies());
  }

private:
  /** Internal type for galaxies. */
  struct internal_galaxy {
    Point point_;
    galaxy_value_type value_;
    /** Constructor. */
    internal_galaxy(const Point& point, const galaxy_value_type& value)
        : point_(point), value_(value) {
    }
  };

  // Store internal galaxy data here.
  std::vector<internal_galaxy> galaxies_;
};
