// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DETAIL_INTERVAL_H_
#define SFCGAL_DETAIL_INTERVAL_H_

#include "SFCGAL/config.h"

namespace SFCGAL::detail {

/**
 * Represents an interval
 */
class SFCGAL_API Interval {
public:
  /**
   * default constructor (empty interval)
   */
  Interval();
  /**
   * collapsed interval constructor
   * @param value The single value for collapsed interval
   */
  Interval(const double &value);
  /**
   * constructor with two values
   * @param value1 First boundary value
   * @param value2 Second boundary value
   */
  Interval(const double &value1, const double &value2);
  /**
   * copy constructor
   * @param other The interval to copy from
   */
  Interval(const Interval &other);
  /**
   * assign operator
   * @param other The interval to assign from
   * @return Reference to this interval
   */
  auto
  operator=(const Interval &other) -> Interval &;

  /**
   * indicates if the interval is empty
   * @return True if the interval is empty
   */
  bool
  isEmpty() const;

  /**
   * return the lower value
   * @return The lower boundary of the interval
   */
  inline const double &
  lower() const
  {
    return _lower;
  }
  /**
   * return the upper value
   * @return The upper boundary of the interval
   */
  inline const double &
  upper() const
  {
    return _upper;
  }
  /**
   * return the width of the interval
   * @return The width (upper - lower) of the interval
   */
  inline double
  width() const
  {
    return _upper - _lower;
  }

  /**
   * expand the interval
   *
   * @warning no effect if isEmpty()
   * @param expandAmount The amount to expand by
   */
  void
  expandBy(const double &expandAmount);
  /**
   * expand the interval to include an other interval.
   *
   * @warning no effect if other.isEmpty()
   * @param other The interval to include
   */
  void
  expandToInclude(const Interval &other);
  /**
   * expand the interval to include a value
   *
   * @warning no effect if value is NaN
   * @param value The value to include
   */
  void
  expandToInclude(const double &value);

  /**
   * test if this intersects other
   * @param other The interval to test intersection with
   * @return True if intervals intersect
   */
  bool
  intersects(const Interval &other) const;

  /**
   * compare two intervals
   * @warning true for empty intervals
   * @param other The interval to compare with
   * @return True if intervals are equal
   */
  bool
  operator==(const Interval &other) const;
  /**
   * compare two intervals
   * @warning false for empty intervals
   * @param other The interval to compare with
   * @return True if intervals are not equal
   */
  bool
  operator!=(const Interval &other) const;

private:
  double _lower;
  double _upper;
};

} // namespace SFCGAL::detail

#endif
