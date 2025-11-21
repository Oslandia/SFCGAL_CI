// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_EQUALITY_H
#define SFCGAL_ALGORITHM_EQUALITY_H

#include "SFCGAL/config.h"
#include <string>

namespace SFCGAL {
class Geometry;

namespace algorithm {

/**
 * Manages how geometry comparison will be handled
 */
class EqualityStrictness {
public:
  // NOLINTBEGIN(performance-enum-size)

  /**
   * @brief Binary enum to define what is allowed when comparing geometries
   *
   * Exclusions:
   * - can not have CheckCoverOrPoint and any of InternalPoint* checks
   * - can not have multiple InternalPoint* checks
   */
  enum Flag {
    /// kind of check: point or cover. 1 is cover, 0 is point
    CheckCoverOrPoint = 1 << 0,
    /// sub geometries need to be ordered
    SubGeomOrdered = 1 << 1,
    /// sub parts need to be ordered
    SubPartOrdered = 1 << 2,
    /// sub geometry points need to be ordered, no shifted is allowed
    InternalPointOrdered = 1 << 3,
    /// sub geometry points need to be ordered but shifted is allowed
    InternalPointShifted = 1 << 4,
    /// sub geometry points need to be ordered but inverted sequence is allowed
    InternalPointInverted = 1 << 5,
  };
  // NOLINTEND(performance-enum-size)

  /// Default constructor
  EqualityStrictness() : _flags(0) {}

  /**
   * Constructor with default value
   * @param flag new flags
   */
  EqualityStrictness(Flag flag) : _flags(flag) {}

  /**
   *  Add flag to active flags
   * @param flag new flag
   * @return ref on this
   */
  auto
  operator|(Flag flag) -> EqualityStrictness &;

  /**
   * Checks if a flag is valid
   * @param flag flag to check
   * @return true when flag is active
   */
  auto
  operator&(Flag flag) const -> bool
  {
    return (_flags & flag) != 0;
  }

  /**
   * Assign operator from Flag value
   * @param flag new flag
   * @return ref on this
   */
  auto
  operator=(Flag flag) -> EqualityStrictness &
  {
    _flags = flag;
    return *this;
  }

  /**
   * Assign operator from EqualityStrictness value
   * @param other another EqualityStrictness object
   * @return ref on this
   */
  auto
  operator=(EqualityStrictness other) -> EqualityStrictness &
  {
    _flags = other._flags;
    return *this;
  }

  /**
   * @return string representation of current flags
   */
  [[nodiscard]] auto
  toString() const -> std::string;

  /**
   * point comparison, sub geometries need to be ordered as sub geometry points
   * @return new object
   */
  static auto
  allPointOrdered() -> EqualityStrictness
  {
    EqualityStrictness out;
    out = out | SubGeomOrdered | SubPartOrdered | InternalPointOrdered;
    return out;
  }

  /**
   * point comparison, nothing need to be ordered
   * @return new object
   */
  static auto
  pointNonOrdered() -> EqualityStrictness
  {
    EqualityStrictness out;
    out._flags = 0;
    return out;
  }

  /**
   * cover comparison, sub geometries do not need to be ordered
   * @return new object
   */
  static auto
  coverSubGeomNonOrdered() -> EqualityStrictness
  {
    EqualityStrictness out;
    out = out | CheckCoverOrPoint;
    return out;
  }

  /**
   * cover comparison, sub geometries need to be ordered
   * @return new object
   */
  static auto
  coverSubGeomOrdered() -> EqualityStrictness
  {
    EqualityStrictness out;
    out = out | CheckCoverOrPoint | SubGeomOrdered;
    return out;
  }

protected:
  /// Holds all flags
  int _flags;
};

/**
 * Equality operator
 * @pre the two geometries must be valid
 * @param geomA geometry to compare with
 * @param geomB geometry to compare with
 * @param tolerance allowed distance between same points.
 * @param strictness will refer to EqualityStrictness flags. Default:
 * EqualityStrictness::allPointOrdered()
 * @return true  when 2 geometries are valid against strictness value and
 * tolerance
 */
SFCGAL_API auto
almostEqual(const Geometry &geomA, const Geometry &geomB, double tolerance,
            EqualityStrictness strictness =
                EqualityStrictness::allPointOrdered()) -> bool;

} // namespace algorithm
} // namespace SFCGAL

#endif // EQUALITY_H
