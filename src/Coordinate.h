// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _SFCGAL_COORDINATE_H_
#define _SFCGAL_COORDINATE_H_

#include <SFCGAL/config.h>

#include <boost/array.hpp>
#include <boost/assert.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/variant.hpp>

#include <SFCGAL/numeric.h>

#include <SFCGAL/Kernel.h>

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
   */
  Coordinate(const Kernel::FT &x, const Kernel::FT &y);
  /**
   * XYZ Constructor with exact coordinates
   */
  Coordinate(const Kernel::FT &x, const Kernel::FT &y, const Kernel::FT &z);
  /**
   * XYZ constructor
   * @warning x,y,z must not be not be NaN nor inf
   */
  Coordinate(const double &x, const double &y, const double &z);

  /**
   * XY constructor
   * @warning x,y must not be not be NaN nor inf
   */
  Coordinate(const double &x, const double &y);
  /**
   * Constructor from CGAL::Point_2<K>
   */
  Coordinate(const Kernel::Point_2 &other);
  /**
   * Constructor from CGAL::Point_3<K>
   */
  Coordinate(const Kernel::Point_3 &other);

  /**
   * copy constructor
   */
  Coordinate(const Coordinate &other);
  /**
   * assign operator
   */
  Coordinate &
  operator=(const Coordinate &other);
  /**
   * destructor
   */
  ~Coordinate();

  /**
   * @brief Get the dimension of the coordinates
   */
  int
  coordinateDimension() const;
  /**
   * @brief Tests if the coordinates are empty
   */
  bool
  isEmpty() const;
  /**
   * @brief Tests if Z is defined
   */
  bool
  is3D() const;

  //--- accessors

  /**
   * @brief Gets the x value
   * @warning Exact, NaN for empty coordinates
   */
  Kernel::FT
  x() const;

  /**
   * @brief Gets the y value
   * @warning Exact, NaN for empty coordinates
   */
  Kernel::FT
  y() const;

  /**
   * @brief Gets the z value
   * @warning Exact, NaN for empty or 0 for 2D coordinates
   */
  Kernel::FT
  z() const;

  //-- helper

  /**
   * @brief round coordinates with a scale factor
   * @return *this
   */
  Coordinate &
  round(const long &scaleFactor = 1);

  //-- comparator

  /**
   * @brief Compares two points (lexicographic order)
   *
   * @warning coordinates must have the same dimension
   */
  bool
  operator<(const Coordinate &other) const;

  /**
   * @brief Compares with an other point
   *
   * @warning coordinates must have the same dimension
   */
  bool
  operator==(const Coordinate &other) const;
  /**
   * @brief Compares with an other point
   *
   * @warning coordinates must have the same dimension
   */
  bool
  operator!=(const Coordinate &other) const;

  /**
   * absolute comparison with an other coordinate
   */
  bool
  almostEqual(const Coordinate &other, const double tolerance) const;

  /**
   * @brief Converts to Kernel::Vector_2
   */
  inline Kernel::Vector_2
  toVector_2() const
  {
    return Kernel::Vector_2(CGAL::ORIGIN, toPoint_2());
  }

  /**
   * @brief Converts to Kernel::Vector_3
   */
  inline Kernel::Vector_3
  toVector_3() const
  {
    return Kernel::Vector_3(CGAL::ORIGIN, toPoint_3());
  }

  /**
   * @brief Converts to Kernel::Point_2
   */
  Kernel::Point_2
  toPoint_2() const;

  /**
   * @brief Converts to Kernel::Point_3
   */
  Kernel::Point_3
  toPoint_3() const;

  // class for Empty coordinate
  class Empty {
  };

private:
  boost::variant<Empty, Kernel::Point_2, Kernel::Point_3> _storage;

public:
  /**
   * Serialization
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

  template <class Archive>
  void
  load(Archive &ar, const unsigned int /*version*/)
  {
    int dim;
    ar >> dim;

    if (dim == 0) {
      _storage = Empty();
    } else if (dim == 2) {
      Kernel::FT x_, y_;
      ar >> x_;
      ar >> y_;
      _storage = Kernel::Point_2(x_, y_);
    } else if (dim == 3) {
      Kernel::FT x_, y_, z_;
      ar >> x_;
      ar >> y_;
      ar >> z_;
      _storage = Kernel::Point_3(x_, y_, z_);
    }
  }

  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version)
  {
    boost::serialization::split_member(ar, *this, version);
  }
};

} // namespace SFCGAL

#endif
