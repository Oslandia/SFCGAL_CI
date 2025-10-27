// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_COORDINATE_H_
#define SFCGAL_COORDINATE_H_

#include "SFCGAL/config.h"

#include <boost/array.hpp>
#include <boost/assert.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/variant.hpp>

#include "SFCGAL/numeric.h"

#include "SFCGAL/Kernel.h"

namespace SFCGAL {

/**
 * @brief Represents the Coordinate of a Point (wraps either an empty structure,
 * or a Kernel::Point_2, or a Kernel::Point_3)
 */
class SFCGAL_API Coordinate {
public:
  /**
   * Empty Coordinate constructor
   */
  Coordinate();

  /**
   * XY Constructor with exact coordinates
   * @param x The X coordinate
   * @param y The Y coordinate
   */
  Coordinate(const Kernel::FT &x, const Kernel::FT &y);
  /**
   * XYZ Constructor with exact coordinates
   * @param x The X coordinate
   * @param y The Y coordinate
   * @param z The Z coordinate
   */
  Coordinate(const Kernel::FT &x, const Kernel::FT &y, const Kernel::FT &z);

  /**
   * XYZ constructor
   * @param x The X coordinate
   * @param y The Y coordinate
   * @param z The Z coordinate
   * @warning x,y,z must not be not be NaN nor inf
   */
  Coordinate(const double &x, const double &y, const double &z);

  /**
   * XY constructor
   * @param x The X coordinate
   * @param y The Y coordinate
   * @warning x,y must not be not be NaN nor inf
   */
  Coordinate(const double &x, const double &y);
  /**
   * Constructor from CGAL::Point_2<K>
   * @param other The CGAL 2D point to copy
   */
  Coordinate(const Kernel::Point_2 &other);
  /**
   * Constructor from CGAL::Point_3<K>
   * @param other The CGAL 3D point to copy
   */
  Coordinate(const Kernel::Point_3 &other);

  /**
   * copy constructor
   * @param other The coordinate to copy from
   */
  Coordinate(const Coordinate &other);
  /**
   * assign operator
   * @param other The coordinate to assign from
   * @return Reference to this coordinate
   */
  auto
  operator=(const Coordinate &other) -> Coordinate &;
  /**
   * destructor
   */
  ~Coordinate();

  /**
   * @brief Get the dimension of the coordinates
   * @return Number of coordinates (2 or 3)
   */
  [[nodiscard]] auto
  coordinateDimension() const -> int;
  /**
   * @brief Tests if the coordinates are empty
   * @return true if empty, false otherwise
   */
  [[nodiscard]] auto
  isEmpty() const -> bool;
  /**
   * @brief Tests if Z is defined
   * @return true if 3D, false otherwise
   */
  [[nodiscard]] auto
  is3D() const -> bool;

  //--- accessors

  /**
   * @brief Gets the x value
   * @return The X coordinate value
   * @warning Exact, NaN for empty coordinates
   */
  [[nodiscard]] auto
  x() const -> Kernel::FT;

  /**
   * @brief Gets the y value
   * @return The Y coordinate value
   * @warning Exact, NaN for empty coordinates
   */
  [[nodiscard]] auto
  y() const -> Kernel::FT;

  /**
   * @brief Gets the z value
   * @return The Z coordinate value
   * @warning Exact, NaN for empty or 0 for 2D coordinates
   */
  [[nodiscard]] auto
  z() const -> Kernel::FT;

  //-- helper

  /**
   * @brief round coordinates with a scale factor
   * @param scaleFactor The scale factor to apply (default 1)
   * @return *this
   */
  auto
  round(const long &scaleFactor = 1) -> Coordinate &;

  //-- comparator

  /**
   * @brief Compares two points (lexicographic order)
   * @param other The coordinate to compare with
   * @return true if this coordinate is less than other
   * @warning coordinates must have the same dimension
   */
  auto
  operator<(const Coordinate &other) const -> bool;

  /**
   * @brief Compares with an other point
   * @param other The coordinate to compare with
   * @return true if coordinates are equal
   * @warning coordinates must have the same dimension
   */
  auto
  operator==(const Coordinate &other) const -> bool;
  /**
   * @brief Compares with an other point
   * @param other The coordinate to compare with
   * @return true if coordinates are not equal
   * @warning coordinates must have the same dimension
   */
  auto
  operator!=(const Coordinate &other) const -> bool;

  /**
   * @brief absolute comparison with an other coordinate
   * @param other The coordinate to compare with
   * @param tolerance The tolerance for comparison
   * @return true if coordinates are almost equal within tolerance
   */
  [[nodiscard]] auto
  almostEqual(const Coordinate &other, double tolerance) const -> bool;

  /**
   * @brief Converts to Kernel::Vector_2
   * @return CGAL 2D vector representation
   */
  [[nodiscard]] auto
  toVector_2() const -> Kernel::Vector_2
  {
    return {CGAL::ORIGIN, toPoint_2()};
  }

  /**
   * @brief Converts to Kernel::Vector_3
   * @return CGAL 3D vector representation
   */
  [[nodiscard]] auto
  toVector_3() const -> Kernel::Vector_3
  {
    return {CGAL::ORIGIN, toPoint_3()};
  }

  /**
   * @brief Converts to Kernel::Point_2
   * @return CGAL 2D point representation
   */
  [[nodiscard]] auto
  toPoint_2() const -> Kernel::Point_2;

  /**
   * @brief Converts to Kernel::Point_3
   * @return CGAL 3D point representation
   */
  [[nodiscard]] auto
  toPoint_3() const -> Kernel::Point_3;

  /**
   * @brief Drops the z coordinate
   * @return TRUE if a Z value was present and has been removed
   */
  auto
  dropZ() -> bool;

  /**
   * @brief Swaps the x and y coordinates
   */
  auto
  swapXY() -> void;

  // class for Empty coordinate
  class Empty {};

private:
  boost::variant<Empty, Kernel::Point_2, Kernel::Point_3> _storage;

public:
  /**
   * @brief Save coordinate to archive
   * @param ar Archive to save to
   */
  template <class Archive>
  void
  save(Archive &ar, const unsigned int /*version*/) const
  {
    int dim = coordinateDimension();
    ar << dim;

    if (_storage.which() > 0) {
      const Kernel::FT &x_ = x();
      const Kernel::FT &y_ = y();
      ar << x_;
      ar << y_;

      if (_storage.which() == 2) {
        const Kernel::FT &z_ = z();
        ar << z_;
      }
    }
  }

  /**
   * @brief Load coordinate from archive
   * @param ar Archive to load from
   */
  template <class Archive>
  void
  load(Archive &ar, const unsigned int /*version*/)
  {
    int dim;
    ar >> dim;

    if (dim == 0) {
      _storage = Empty();
    } else if (dim == 2) {
      Kernel::FT x_;
      Kernel::FT y_;
      ar >> x_;
      ar >> y_;
      _storage = Kernel::Point_2(x_, y_);
    } else if (dim == 3) {
      Kernel::FT x_;
      Kernel::FT y_;
      Kernel::FT z_;
      ar >> x_;
      ar >> y_;
      ar >> z_;
      _storage = Kernel::Point_3(x_, y_, z_);
    }
  }

  /**
   * @brief Serialize coordinate to/from archive
   * @param ar Archive for serialization
   * @param version Serialization version
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version)
  {
    boost::serialization::split_member(ar, *this, version);
  }
};

} // namespace SFCGAL

#endif
