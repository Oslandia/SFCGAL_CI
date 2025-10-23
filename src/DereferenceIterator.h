// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DEREFERENCEITERATOR_H_
#define SFCGAL_DEREFERENCEITERATOR_H_

#include <cstddef>
#include <type_traits>

namespace SFCGAL {

/**
 * @brief Iterator adapter that dereferences unique_ptr elements automatically.
 *
 * This class wraps a base iterator (e.g.,
 * std::vector<std::unique_ptr<T>>::iterator) and returns references to the
 * pointed objects instead of unique_ptrs.
 *
 * @tparam BaseIterator the underlying iterator type
 */
template <class BaseIterator>
class DereferenceIterator : public BaseIterator {
public:
  using value_type =
      typename BaseIterator::value_type::element_type; ///< Dereferenced
  using pointer   = value_type *; ///< Pointer to value type
  using reference = value_type &; ///< Reference to value type

  /** @brief Default constructor. */
  DereferenceIterator() = default;

  /**
   * @brief Construct from a base iterator.
   * @param other The underlying base iterator.
   */
  explicit DereferenceIterator(BaseIterator other) : BaseIterator(other) {}

  /**
   * @brief Conversion constructor from another compatible DereferenceIterator.
   *
   * Allows implicit conversion from DereferenceIterator<iterator>
   * to DereferenceIterator<const_iterator>.
   *
   * @tparam OtherIterator Another iterator type convertible to BaseIterator.
   * @param other The other DereferenceIterator to copy.
   */
  template <typename OtherIterator,
            typename = std::enable_if_t<
                std::is_convertible_v<OtherIterator, BaseIterator>>>
  DereferenceIterator(const DereferenceIterator<OtherIterator> &other)
      : BaseIterator(other)
  {
  }

  /**
   * @brief Dereference operator returning a reference to the pointed value.
   * @return Reference to the object pointed by the current iterator.
   */
  auto
  operator*() const -> reference
  {
    return *(this->BaseIterator::operator*());
  }

  /**
   * @brief Member access operator returning a pointer to the pointed value.
   * @return Pointer to the object pointed by the current iterator.
   */
  auto
  operator->() const -> pointer
  {
    return this->BaseIterator::operator*().get();
  }

  /**
   * @brief Random access operator returning a reference to the nth value.
   * @param idx The offset from the current iterator position.
   * @return Reference to the nth object from the current iterator.
   */
  auto
  operator[](size_t idx) const -> reference
  {
    return *(this->BaseIterator::operator[](idx));
  }
};

/**
 * @brief Helper function to create a DereferenceIterator from a base iterator.
 *
 * @tparam Iterator The underlying iterator type.
 * @param iterator The base iterator to wrap.
 * @return A DereferenceIterator instance wrapping @p t.
 */
template <typename Iterator>
auto
dereference_iterator(Iterator iterator) -> DereferenceIterator<Iterator>
{
  return DereferenceIterator<Iterator>(iterator);
}

} // namespace SFCGAL

#endif // SFCGAL_DEREFERENCEITERATOR_H_
