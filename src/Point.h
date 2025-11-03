// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_POINT_H_
#define SFCGAL_POINT_H_

#include "SFCGAL/Coordinate.h"

#include "SFCGAL/Geometry.h"
#include "SFCGAL/detail/TypeForDimension.h"

#include <boost/serialization/base_object.hpp>

namespace SFCGAL {

/**
 * A point in SFA.
 * The x(),y(),z() interface is based on CGAL kernel requirements, taken
 * from examples/Kernel_23/MyPointC2.h
 * @todo strong typing on coordinate dimension?
 */
class SFCGAL_API Point : public GeometryImpl<Point, Geometry> {
public:
  /**
   * Empty point constructor
   */
  Point();
  /**
   * Constructor with Coordinate
   * @param coordinate The coordinate to use for this point
   */
  Point(const Coordinate &coordinate);
  /**
   * XY Constructor with exact coordinates
   * @param x The X coordinate
   * @param y The Y coordinate
   */
  Point(const Kernel::FT &x, const Kernel::FT &y);
  /**
   * XYZ Constructor with exact coordinates
   * @param x The X coordinate
   * @param y The Y coordinate
   * @param z The Z coordinate
   * @param m The M coordinate (default NaN)
   */
  Point(const Kernel::FT &x, const Kernel::FT &y, const Kernel::FT &z,
        const double &m = NaN());
  /**
   * XY constructor
   * @param x The X coordinate
   * @param y The Y coordinate
   */
  Point(const double &x, const double &y);

  /**
   * XYZ constructor
   * @param x The X coordinate
   * @param y The Y coordinate
   * @param z The Z coordinate
   */
  Point(const double &x, const double &y, const double &z);

  /**
   * XYZM constructor
   * @param x The X coordinate
   * @param y The Y coordinate
   * @param z The Z coordinate
   * @param m The M coordinate
   */
  Point(const double &x, const double &y, const double &z, const double &m);

  /**
   * Generic XYZM constructor determined by CoordinateType
   * @param x The X coordinate
   * @param y The Y coordinate
   * @param z The Z coordinate
   * @param m The M coordinate
   * @param dim The coordinate type dimension
   */
  Point(const double &x, const double &y, const double &z, const double &m,
        CoordinateType dim);

  /**
   * Generic XYZM constructor determined by CoordinateType
   * @param x The X coordinate
   * @param y The Y coordinate
   * @param z The Z coordinate
   * @param m The M coordinate
   * @param dim The coordinate type dimension
   */
  Point(const Kernel::FT &x, const Kernel::FT &y, const Kernel::FT &z,
        const double &m, CoordinateType dim);

  /**
   * Constructor from CGAL::Point_2<K>
   * @param other The CGAL 2D point to copy
   */
  Point(const Kernel::Point_2 &other);
  /**
   * Constructor from CGAL::Point_3<K>
   * @param other The CGAL 3D point to copy
   */
  Point(const Kernel::Point_3 &other);

  /**
   * copy constructor
   * @param other The point to copy from
   */
  Point(const Point &other);
  /**
   * assign operator
   * @param other The point to assign from
   * @return Reference to this point
   */
  auto
  operator=(const Point &other) -> Point &;
  /**
   * destructor
   */
  ~Point() override;

  //-- SFCGAL::Geometry
  /// @brief Get the geometry type as string
  /// @return "Point"
  [[nodiscard]] auto
  geometryType() const -> std::string override;
  //-- SFCGAL::Geometry
  /// @brief Get the geometry type identifier
  /// @return TYPE_POINT
  [[nodiscard]] auto
  geometryTypeId() const -> GeometryType override;
  //-- SFCGAL::Geometry
  /// @brief Get the dimension of the point
  /// @return 0 (points are 0-dimensional)
  [[nodiscard]] auto
  dimension() const -> int override;
  //-- SFCGAL::Geometry
  /// @brief Get the coordinate dimension
  /// @return Number of coordinates (2, 3, or 4)
  [[nodiscard]] auto
  coordinateDimension() const -> int override;
  //-- SFCGAL::Geometry
  /// @brief Check if the point is empty
  /// @return true if empty, false otherwise
  [[nodiscard]] auto
  isEmpty() const -> bool override;
  //-- SFCGAL::Geometry
  /// @brief Check if the point has 3D coordinates
  /// @return true if 3D, false otherwise
  [[nodiscard]] auto
  is3D() const -> bool override;
  //-- SFCGAL::Geometry
  /// @brief Check if the point has measured coordinates
  /// @return true if measured, false otherwise
  [[nodiscard]] auto
  isMeasured() const -> bool override;

  /// @brief Drop Z coordinate
  /// @return true if Z was dropped, false otherwise
  auto
  dropZ() -> bool override;

  /// @brief Drop M coordinate
  /// @return true if M was dropped, false otherwise
  auto
  dropM() -> bool override;

  /// @brief Swap X and Y coordinates
  auto
  swapXY() -> void override;

  //--- accessors

  /**
   * @brief Returns the x value as a double throw for empty Point
   * @return X coordinate value
   */
  [[nodiscard]] auto
  x() const -> Kernel::RT
  {
    return _coordinate.x();
  }
  /**
   * @brief Returns the y value as a double throw for empty Point
   * @return Y coordinate value
   */
  [[nodiscard]] auto
  y() const -> Kernel::RT
  {
    return _coordinate.y();
  }
  /**
   * @brief Returns the z value (zero for 2D) throw for empty Point
   * @return Z coordinate value
   */
  [[nodiscard]] auto
  z() const -> Kernel::RT
  {
    return _coordinate.z();
  }

  /**
   * @brief Returns the m value (NaN is not defined)
   * @return M coordinate value
   */
  [[nodiscard]] auto
  m() const -> double
  {
    return _m;
  }
  /**
   * Sets the m value
   * @param m The new M coordinate value
   */
  void
  setM(const double &m)
  {
    _m = m;
  }

  /**
   * compare two points
   * @param other The point to compare with
   * @return true if this point is less than other
   */
  auto
  operator<(const Point &other) const -> bool;

  /**
   * compare with an other point
   * @param other The point to compare with
   * @return true if points are equal
   */
  auto
  operator==(const Point &other) const -> bool;
  /**
   * compare with an other point
   * @param other The point to compare with
   * @return true if points are not equal
   */
  auto
  operator!=(const Point &other) const -> bool;

  /**
   * absolute comparison with an other point
   * @param other The point to compare with
   * @param tolerance The tolerance for comparison
   * @return true if points are almost equal within tolerance
   */
  [[nodiscard]] auto
  almostEqual(const Point &other, double tolerance) const -> bool;

  //-- visitors

  //-- SFCGAL::Geometry
  /// @brief Accept a geometry visitor
  /// @param visitor Visitor to accept
  void
  accept(GeometryVisitor &visitor) override;
  //-- SFCGAL::Geometry
  /// @brief Accept a const geometry visitor
  /// @param visitor Const visitor to accept
  void
  accept(ConstGeometryVisitor &visitor) const override;

  /**
   * @brief Convert to CGAL 2D vector
   * @return CGAL 2D vector representation
   * @see Coordinate::toVector_2()
   */
  [[nodiscard]] auto
  toVector_2() const -> Kernel::Vector_2
  {
    return _coordinate.toVector_2();
  }

  /**
   * @brief Convert to CGAL 3D vector
   * @return CGAL 3D vector representation
   * @see Coordinate::toVector_3()
   */
  [[nodiscard]] auto
  toVector_3() const -> Kernel::Vector_3
  {
    return _coordinate.toVector_3();
  }

  /**
   * @brief Convert to CGAL 2D point
   * @return CGAL 2D point representation
   * @see Coordinate::toPoint_2()
   */
  [[nodiscard]] auto
  toPoint_2() const -> Kernel::Point_2
  {
    return _coordinate.toPoint_2();
  }

  /**
   * @brief Convert to CGAL 3D point
   * @return CGAL 3D point representation
   * @see Coordinate::toPoint_3()
   */
  [[nodiscard]] auto
  toPoint_3() const -> Kernel::Point_3
  {
    return _coordinate.toPoint_3();
  }

  /**
   * @brief Converts to CGAL::Point_2 or CGAL::Point_3
   * @tparam Dim The dimension (2 or 3)
   * @return CGAL point of specified dimension
   */
  template <int Dim>
  auto
  toPoint_d() const -> typename detail::TypeForDimension<Dim>::Point;

  /**
   * @brief Get reference to coordinate
   * @return Reference to the coordinate
   */
  auto
  coordinate() -> Coordinate &
  {
    return _coordinate;
  }
  /**
   * @brief Get const reference to coordinate
   * @return Const reference to the coordinate
   */
  [[nodiscard]] auto
  coordinate() const -> const Coordinate &
  {
    return _coordinate;
  }

  /**
   * Serializer
   * @param ar Archive for serialization
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int /*version*/)
  {
    ar &boost::serialization::base_object<Geometry>(*this);
    ar & _coordinate;
  }

private:
  template <typename T>
  void
             newPoint(const T &x, const T &y, const T &z, const double &m,
                      CoordinateType dim);
  Coordinate _coordinate;
  /**
   * @brief m coordinates (NaN if not defined)
   */
  double _m;
};

} // namespace SFCGAL

#endif
