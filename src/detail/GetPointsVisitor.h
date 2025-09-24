// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DETAIL_GETPOINTSVISITOR_H_
#define SFCGAL_DETAIL_GETPOINTSVISITOR_H_

#include "SFCGAL/config.h"

#include "SFCGAL/GeometryVisitor.h"
#include <vector>

namespace SFCGAL {
namespace detail {

/**
 * Get the list of points from a Geometry
 */
class SFCGAL_API GetPointsVisitor : public ConstGeometryVisitor {
public:
  /**
   * @brief Visit a Point and add it to collection
   * @param g The point to visit
   */
  void
  visit(const Point &g) override;
  /**
   * @brief Visit a LineString and add all points to collection
   * @param g The linestring to visit
   */
  void
  visit(const LineString &g) override;
  /**
   * @brief Visit a Polygon and add all points to collection
   * @param g The polygon to visit
   */
  void
  visit(const Polygon &g) override;
  /**
   * @brief Visit a Triangle and add all points to collection
   * @param g The triangle to visit
   */
  void
  visit(const Triangle &g) override;
  /**
   * @brief Visit a Solid and add all points to collection
   * @param g The solid to visit
   */
  void
  visit(const Solid &g) override;
  /**
   * @brief Visit a MultiPoint and add all points to collection
   * @param g The multipoint to visit
   */
  void
  visit(const MultiPoint &g) override;
  /**
   * @brief Visit a MultiLineString and add all points to collection
   * @param g The multilinestring to visit
   */
  void
  visit(const MultiLineString &g) override;
  /**
   * @brief Visit a MultiPolygon and add all points to collection
   * @param g The multipolygon to visit
   */
  void
  visit(const MultiPolygon &g) override;
  /**
   * @brief Visit a MultiSolid and add all points to collection
   * @param g The multisolid to visit
   */
  void
  visit(const MultiSolid &g) override;
  /**
   * @brief Visit a GeometryCollection and add all points to collection
   * @param g The geometry collection to visit
   */
  void
  visit(const GeometryCollection &g) override;
  /**
   * @brief Visit a PolyhedralSurface and add all points to collection
   * @param g The polyhedral surface to visit
   */
  void
  visit(const PolyhedralSurface &g) override;
  /**
   * @brief Visit a TriangulatedSurface and add all points to collection
   * @param g The triangulated surface to visit
   */
  void
  visit(const TriangulatedSurface &g) override;

public:
  typedef std::vector<const Point *>::const_iterator
                             const_iterator; ///< Const iterator type for points
  std::vector<const Point *> points;         ///< Collection of collected points
};

} // namespace detail
} // namespace SFCGAL

#endif
