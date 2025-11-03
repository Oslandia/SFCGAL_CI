// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DETAIL_FORCEVALIDITY_VISITOR_H_
#define SFCGAL_DETAIL_FORCEVALIDITY_VISITOR_H_

#include "SFCGAL/config.h"

#include "SFCGAL/GeometryVisitor.h"

namespace SFCGAL::detail {

/**
 * @brief Visitor that forces the validity flag on geometries
 *
 * This visitor traverses a geometry hierarchy and sets the validity
 * flag to the specified value for all visited geometries.
 */
class SFCGAL_API ForceValidityVisitor : public GeometryVisitor {
public:
  /**
   * @brief Construct validity visitor
   * @param valid The validity flag to apply to all visited geometries
   */
  ForceValidityVisitor(bool valid) : valid_(valid) {}
  /**
   * @brief Visit and set validity flag on Point
   * @param g The point to visit
   */
  void
  visit(Point &g) override;
  /**
   * @brief Visit and set validity flag on LineString
   * @param g The linestring to visit
   */
  void
  visit(LineString &g) override;
  /**
   * @brief Visit and set validity flag on Polygon
   * @param g The polygon to visit
   */
  void
  visit(Polygon &g) override;
  /**
   * @brief Visit and set validity flag on Triangle
   * @param g The triangle to visit
   */
  void
  visit(Triangle &g) override;
  /**
   * @brief Visit and set validity flag on Solid
   * @param g The solid to visit
   */
  void
  visit(Solid &g) override;
  /**
   * @brief Visit and set validity flag on MultiPoint
   * @param g The multipoint to visit
   */
  void
  visit(MultiPoint &g) override;
  /**
   * @brief Visit and set validity flag on MultiLineString
   * @param g The multilinestring to visit
   */
  void
  visit(MultiLineString &g) override;
  /**
   * @brief Visit and set validity flag on MultiPolygon
   * @param g The multipolygon to visit
   */
  void
  visit(MultiPolygon &g) override;
  /**
   * @brief Visit and set validity flag on MultiSolid
   * @param g The multisolid to visit
   */
  void
  visit(MultiSolid &g) override;
  /**
   * @brief Visit and set validity flag on GeometryCollection
   * @param g The geometry collection to visit
   */
  void
  visit(GeometryCollection &g) override;
  /**
   * @brief Visit and set validity flag on PolyhedralSurface
   * @param g The polyhedral surface to visit
   */
  void
  visit(PolyhedralSurface &g) override;
  /**
   * @brief Visit and set validity flag on TriangulatedSurface
   * @param g The triangulated surface to visit
   */
  void
  visit(TriangulatedSurface &g) override;
  /**
   * @brief Visit and set validity flag on NURBSCurve
   * @param g The triangulated surface to visit
   */
  void
  visit(NURBSCurve &g) override;

private:
  bool valid_;
};

} // namespace SFCGAL::detail

#endif
