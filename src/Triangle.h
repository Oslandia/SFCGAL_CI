// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TRIANGLE_H_
#define SFCGAL_TRIANGLE_H_

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
class SFCGAL_API Triangle : public Surface {
public:
  /**
   * empty Triangle constructor
   */
  Triangle();
  /**
   * Constructor with a CGAL triangle
   */
  Triangle(const Kernel::Triangle_2 &triangle);
  /**
   * Constructor with a CGAL triangle
   */
  Triangle(const Kernel::Triangle_3 &triangle);
  /**
   * constructor with 3 points
   */
  Triangle(const Point &p, const Point &q, const Point &r);
  /**
   * copy constructor
   */
  Triangle(const Triangle &other);
  /**
   * assign operator
   */
  Triangle &
  operator=(const Triangle &other);
  /**
   * destructor
   */
  ~Triangle();

  //-- SFCGAL::Geometry
  virtual Triangle *
  clone() const;

  //-- SFCGAL::Geometry
  virtual std::string
  geometryType() const;
  //-- SFCGAL::Geometry
  virtual GeometryType
  geometryTypeId() const;
  //-- SFCGAL::Geometry
  virtual int
  coordinateDimension() const;
  //-- SFCGAL::Geometry
  virtual bool
  isEmpty() const;
  //-- SFCGAL::Geometry
  virtual bool
  is3D() const;
  //-- SFCGAL::Geometry
  virtual bool
  isMeasured() const;

  /**
   * reverse Triangle orientation
   */
  void
  reverse();

  /**
   * convert a triangle to a polygon
   */
  Polygon
  toPolygon() const;

  /**
   * returns the i-th vertex
   */
  inline const Point &
  vertex(const int &i) const
  {
    return _vertices[i % 3];
  }
  /**
   * returns the i-th vertex
   */
  inline Point &
  vertex(const int &i)
  {
    return _vertices[i % 3];
  }

  /**
   * Convert to CGAL::Triangle_2
   */
  inline Kernel::Triangle_2
  toTriangle_2() const
  {
    return Kernel::Triangle_2(vertex(0).toPoint_2(), vertex(1).toPoint_2(),
                              vertex(2).toPoint_2());
  }

  /**
   * Convert to CGAL::Triangle_3
   */
  inline Kernel::Triangle_3
  toTriangle_3() const
  {
    return Kernel::Triangle_3(vertex(0).toPoint_3(), vertex(1).toPoint_3(),
                              vertex(2).toPoint_3());
  }

  /**
   * Convert to CGAL::Triangle_2 or CGAL::Triangle_2
   */
  template <int D>
  inline typename detail::TypeForDimension<D>::Triangle
  toTriangle_d() const
  {
    return typename detail::TypeForDimension<D>::Triangle(
        vertex(0).toPoint_d<D>(), vertex(1).toPoint_d<D>(),
        vertex(2).toPoint_d<D>());
  }

  //-- visitors

  //-- SFCGAL::Geometry
  virtual void
  accept(GeometryVisitor &visitor);
  //-- SFCGAL::Geometry
  virtual void
  accept(ConstGeometryVisitor &visitor) const;

  /**
   * Serializer
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
  Point _vertices[3];
};

} // namespace SFCGAL

#endif
