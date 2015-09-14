#pragma once
/** @file SpaceSearcher.hpp
 * @brief Define the SpaceSearcher class for making spatial searches.
 */

#include "Point.hpp"
#include "BoundingBox.hpp"
#include "MortonCoder.hpp"
#include "../Util/GeneralUtil.hpp"  // totally_ordered

/** @class SpaceSearcher
 * @brief Class for making spatial searches, which uses the MortonCoder
 *        class as a backend.
 *
 * Given a range of data items (e.g., Nodes) and a mapping between these
 * data items and Points, the SpaceSearcher class can be used to quickly
 * iterate over data items which are contained within an epsilon of
 * any given BoundingBox.
 *
 * See "morton_test.cpp" for an usage example.
 */
template <typename T, typename T2Point, int L = 10>
class SpaceSearcher {
 public:

  ////////////////////////////////////
  // TYPE DEFINITIONS AND CONSTANTS //
  ////////////////////////////////////

  /** The number of levels in the MortonCoder. */
  static constexpr int NumLevels = L;
  /** Type of MortonCoder. */
  typedef MortonCoder<NumLevels> MortonCoderType;
  /** Type of the Morton codes. */
  typedef typename MortonCoderType::code_type code_type;

  /** Type of iterators, which iterate over items inside a BoundingBox. */
  class Iterator;

  /** Synonym for Iterator */
  typedef Iterator iterator;
  typedef Iterator const_iterator;

 private:
  // Implementation types
  using pair = std::pair<code_type,T>;

 public:

  /////////////////
  // CONSTRUCTOR //
  /////////////////

  /** @brief Constructor.
   *
   * For a range of data items of type @a T given by [@a t_begin, @a t_end)
   * and a mapping @a t2p between data items and Points, initialize a
   * framework for making spatial searches.
   *
   * @param[in] t_begin Iterator to first data item.
   * @param[in] t_end   Iterator to one past the last data item.
   * @param[in] t2p     A functor that maps data items to Points.
   */
  template <typename TIter>
  SpaceSearcher(TIter t_begin, TIter t_end, T2Point t2p) {
    // Static if reserve optimization
    if (std::is_same<typename std::iterator_traits<TIter>::iterator_category,
                     std::random_access_iterator_tag>::value)
      t_.reserve(std::distance(t_begin,t_end));

    // Construct the enclosing bounding box and copy the T
    BoundingBox bb;
    for (auto it = t_begin; it != t_end; ++it) {
      t_.emplace_back(0,*it);
      bb |= t2p(t_.back().second);
    }

    // Create MortonCoder instance.
    mc_ = MortonCoderType(bb);

    // Fill in the morton codes
    for (auto& p : t_)
      p.first = mc_.code(t2p(p.second));

    // Sort by morton code
    std::sort(t_.begin(), t_.end(),
              [](const pair& a, const pair& b) { return a.first < b.first; });
  }

  //////////////
  // Iterator //
  //////////////

  /** @class SpaceSearcher::Iterator
   * @brief Iterator class for data items. A forward iterator.
   *
   * Iterates over data items of type @a T contained
   * within epsilon of a given bounding box.
   */
  class Iterator {
   public:
    // These type definitions help us use STL's iterator_traits.
    typedef T value_type;
    typedef T* pointer;
    typedef T& reference;
    typedef std::ptrdiff_t difference_type;
    typedef std::forward_iterator_tag iterator_category;

    // Default constructor
    Iterator() = default;

    const T& operator*() const {
      return i_->second;
    }
    Iterator& operator++() {
      ++i_; fix();
      return *this;
    }
    bool operator==(const Iterator& other) const {
      return i_ == other.i_;
    }
    bool operator!=(const Iterator& other) const {
      return i_ != other.i_;
    }

   private:
    friend SpaceSearcher;
    using m_iter = typename std::vector<pair>::const_iterator;
    m_iter i_, end_;
    code_type min_, max_;
    Iterator(m_iter i, m_iter end, code_type min, code_type max)
        : i_(i), end_(end), min_(min), max_(max) {
      fix();
    }
    // @post mi_ == end_ || MortonCoderType::is_in_box(mi_->first, min_, max_)
    void fix() {
      while (i_ < end_) {
        code_type c = MortonCoderType::advance_to_box(i_->first, min_, max_);
        if (i_->first == c)
          break;
        i_ = std::partition_point(i_, end_, [&](const pair& p) {
            return p.first < c;
          });
      }
    }
  };

  /** Iterator to the first item
   * contained within some epsilon of a bounding box.
   */
  const_iterator begin(BoundingBox bb) const {
    // Intersection of bb and our domain
    bb &= mc_.bounding_box();
    code_type morton_min = mc_.code(bb.min());
    code_type morton_max = mc_.code(bb.max());
    auto mit = std::partition_point(t_.begin(), t_.end(),
                                    [&](const pair& p) {
                                      return p.first < morton_min;
                                    });
    auto mit_end = std::partition_point(mit, t_.end(),
                                        [&](const pair& p) {
                                          return p.first < morton_max;
                                        });
    return Iterator(mit, mit_end, morton_min, morton_max);
  }

  /** Iterator to one-past-the-last item
   * contained within some epsilon of a bounding box
   */
  const_iterator end(BoundingBox bb) const {
    // Intersection of bb and our domain
    bb &= mc_.bounding_box();
    //code_type morton_min = mc_.code(bb.min());
    code_type morton_max = mc_.code(bb.max());
    auto mit_end = std::partition_point(t_.begin(), t_.end(),
                                        [&](const pair& p) {
                                          return p.first < morton_max;
                                        });
    return Iterator(mit_end, mit_end, code_type(0), code_type(0));
  }

 private:
  // MortonCoder instance associated with this SpaceSearcher.
  MortonCoderType mc_;
  // Pairs of Morton codes and data items of type T.
  std::vector<std::pair<code_type,T>> t_;
};
