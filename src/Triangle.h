// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TRIANGLE_H_
#define SFCGAL_TRIANGLE_H_

#include <array>
#include <boost/shared_ptr.hpp>

#include <boost/serialization/base_object.hpp>

#include "SFCGAL/Point.h"
#include "SFCGAL/Surface.h"
#include "SFCGAL/detail/TypeForDimension.h"

#include <CGAL/Triangle_2.h>
#include <CGAL/Triangle_3.h>

namespace SFCGAL {

/**
 * [OGC/SFA]Triangle
 *
 * @warning According to SFA, a Triangle should be inherited from a Polygon.
 * That means that a triangle "is a" Polygon with hole. This inheritance is
 * removed in order to keep CGAL modeling.
 *
 * @warning An empty triangle has empty points
 */
class SFCGAL_API Triangle : public GeometryImpl<Triangle, Surface> {
public:
  /**
   * empty Triangle constructor
   */
  Triangle();
  /**
   * @brief Constructor with a CGAL triangle
   * @param triangle CGAL 2D triangle to construct from
   */
  Triangle(const Kernel::Triangle_2 &triangle);
  /**
   * @brief Constructor with a CGAL triangle
   * @param triangle CGAL 3D triangle to construct from
   */
  Triangle(const Kernel::Triangle_3 &triangle);
  /**
   * constructor with 3 points
   * @param point1 First point of the triangle
   * @param point2 Second point of the triangle
   * @param point3 Third point of the triangle
   */
  Triangle(const Point &point1, const Point &point2, const Point &point3);
  /**
   * copy constructor
   * @param other Triangle to copy from
   */
  Triangle(const Triangle &other);
  /**
   * @brief assign operator
   * @param other Triangle to assign from
   * @return Reference to this triangle
   */
  auto
  operator=(const Triangle &other) -> Triangle &;
  /**
   * destructor
   */
  ~Triangle() override;

  //-- SFCGAL::Geometry
  /// @brief Get the geometry type as string
  /// @return Geometry type string
  [[nodiscard]] auto
  geometryType() const -> std::string override;
  //-- SFCGAL::Geometry
  /// @brief Get the geometry type identifier
  /// @return Geometry type ID
  [[nodiscard]] auto
  geometryTypeId() const -> GeometryType override;
  //-- SFCGAL::Geometry
  /// @brief Get the coordinate dimension
  /// @return Coordinate dimension
  [[nodiscard]] auto
  coordinateDimension() const -> int override;
  //-- SFCGAL::Geometry
  /// @brief Check if the triangle is empty
  /// @return true if empty, false otherwise
  [[nodiscard]] auto
  isEmpty() const -> bool override;
  //-- SFCGAL::Geometry
  /// @brief Check if the triangle is 3D
  /// @return true if 3D, false otherwise
  [[nodiscard]] auto
  is3D() const -> bool override;
  //-- SFCGAL::Geometry
  /// @brief Check if the triangle has measured coordinates
  /// @return true if measured, false otherwise
  [[nodiscard]] auto
  isMeasured() const -> bool override;

  /// @brief Drop Z coordinate from all points
  /// @return true if Z was dropped, false otherwise
  auto
  dropZ() -> bool override;

  /// @brief Drop M coordinate from all points
  /// @return true if M was dropped, false otherwise
  auto
  dropM() -> bool override;

  /// @brief Swap X and Y coordinates
  auto
  swapXY() -> void override;

  /**
   * reverse Triangle orientation
   */
  void
  reverse();

  /**
   * @brief convert a triangle to a polygon
   * @return Polygon representation of the triangle
   */
  [[nodiscard]] auto
  toPolygon() const -> Polygon;

  /**
   * returns the i-th vertex
   * @param i The vertex index (modulo 3)
   * @return Const reference to the vertex
   */
  [[nodiscard]] auto
  vertex(const int &i) const -> const Point &
  {
    return _vertices[i % 3];
  }
  /**
   * returns the i-th vertex
   * @param i The vertex index (modulo 3)
   * @return Reference to the vertex
   */
  auto
  vertex(const int &i) -> Point &
  {
    return _vertices[i % 3];
  }

  /**
   * Convert to CGAL::Triangle_2
   * @return 2D triangle representation
   */
  [[nodiscard]] auto
  toTriangle_2() const -> Kernel::Triangle_2
  {
    return {vertex(0).toPoint_2(), vertex(1).toPoint_2(),
            vertex(2).toPoint_2()};
  }

  /**
   * Convert to CGAL::Triangle_3
   * @return 3D triangle representation
   */
  [[nodiscard]] auto
  toTriangle_3() const -> Kernel::Triangle_3
  {
    return {vertex(0).toPoint_3(), vertex(1).toPoint_3(),
            vertex(2).toPoint_3()};
  }

  /**
   * Convert to CGAL::Triangle_2 or CGAL::Triangle_3
   * @tparam D Dimension (2 or 3)
   * @return Triangle in specified dimension
   */
  template <int D>
  auto
  toTriangle_d() const -> typename detail::TypeForDimension<D>::Triangle
  {
    return typename detail::TypeForDimension<D>::Triangle(
        vertex(0).toPoint_d<D>(), vertex(1).toPoint_d<D>(),
        vertex(2).toPoint_d<D>());
  }

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
   * Serializer
   * @param ar Archive for serialization
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int /*version*/)
  {
    ar &boost::serialization::base_object<Geometry>(*this);
    ar &_vertices[0] & _vertices[1] & _vertices[2];
  }

private:
  /**
   * point forming the triangle
   */
  std::array<Point, 3> _vertices;
};

} // namespace SFCGAL

#endif
