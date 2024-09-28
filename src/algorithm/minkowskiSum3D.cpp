// Copyright (c) 2012-2024, SFCGAL Contributors and Oslandia
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/minkowskiSum3D.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/detail/GeometrySet.h"
#include "SFCGAL/detail/transform/AffineTransform3.h"

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/minkowski_sum_3.h>

namespace SFCGAL::algorithm {

using Nef_polyhedron_3 = CGAL::Nef_polyhedron_3<Kernel>;
using Polyhedron_3     = CGAL::Polyhedron_3<Kernel>;

auto
perpendicular_vector(const Kernel::Vector_3 &v) -> Kernel::Vector_3
{
  if (v.x() != 0 || v.y() != 0) {
    return Kernel::Vector_3(-v.y(), v.x(), 0);
  }     return Kernel::Vector_3(0, -v.z(), v.y());
 
}

// Helper function to convert SFCGAL::Geometry to Nef_polyhedron_3
auto
geometryToNef(const Geometry &g) -> Nef_polyhedron_3
{
  Nef_polyhedron_3                      result;
  std::function<void(const Geometry &)> process_geometry =
      [&](const Geometry &geom) {
        switch (geom.geometryTypeId()) {
        case TYPE_POINT: {
          const auto &p = geom.as<Point>();
          result        = Nef_polyhedron_3(p.toPoint_3());
          break;
        }
        case TYPE_LINESTRING: {
          const auto &ls = geom.as<LineString>();
          if (ls.numPoints() < 2) {
            result = Nef_polyhedron_3();
          } else {
            std::vector<Kernel::Point_3> points;
            for (size_t i = 0; i < ls.numPoints(); ++i) {
              points.push_back(ls.pointN(i).toPoint_3());
            }
            Kernel::FT   radius(0.001);
            Polyhedron_3 poly;
            for (size_t i = 0; i < points.size() - 1; ++i) {
              Kernel::Vector_3 dir   = points[i + 1] - points[i];
              Kernel::Vector_3 perp1 = perpendicular_vector(dir);
              Kernel::Vector_3 perp2 = CGAL::cross_product(dir, perp1);

              double perp1_length =
                  std::sqrt(CGAL::to_double(perp1.squared_length()));
              double perp2_length =
                  std::sqrt(CGAL::to_double(perp2.squared_length()));

              perp1 = perp1 * (radius / Kernel::FT(perp1_length));
              perp2 = perp2 * (radius / Kernel::FT(perp2_length));

              std::vector<Kernel::Point_3> segment_points = {
                  points[i] + perp1,     points[i] - perp1,
                  points[i] + perp2,     points[i] - perp2,
                  points[i + 1] + perp1, points[i + 1] - perp1,
                  points[i + 1] + perp2, points[i + 1] - perp2};
              CGAL::convex_hull_3(segment_points.begin(), segment_points.end(),
                                  poly);
            }
            result = Nef_polyhedron_3(poly);
          }
          break;
        }
        case TYPE_TRIANGLE: {
          const auto  &tri = geom.as<Triangle>();
          Polyhedron_3 poly;
          poly.make_triangle(tri.vertex(0).toPoint_3(),
                             tri.vertex(1).toPoint_3(),
                             tri.vertex(2).toPoint_3());
          result = Nef_polyhedron_3(poly);
          break;
        }
        case TYPE_POLYGON: {
          const auto                  &poly = geom.as<Polygon>();
          std::vector<Kernel::Point_3> points;
          for (size_t i = 0; i < poly.exteriorRing().numPoints() - 1; ++i) {
            points.push_back(poly.exteriorRing().pointN(i).toPoint_3());
          }
          Polyhedron_3 cgal_poly;
          CGAL::convex_hull_3(points.begin(), points.end(), cgal_poly);
          result = Nef_polyhedron_3(cgal_poly);
          break;
        }
        case TYPE_TRIANGULATEDSURFACE:
        case TYPE_POLYHEDRALSURFACE: {
          const auto  &ps = geom.as<PolyhedralSurface>();
          Polyhedron_3 poly;
          for (size_t i = 0; i < ps.numPolygons(); ++i) {
            Nef_polyhedron_3 temp_nef = geometryToNef(ps.polygonN(i));
            if (i == 0) {
              result = temp_nef;
            } else {
              result += temp_nef;
            }
          }
          break;
        }
        case TYPE_SOLID: {
          const auto &solid = geom.as<Solid>();
          process_geometry(solid.exteriorShell());
          break;
        }
        case TYPE_MULTIPOINT:
        case TYPE_MULTILINESTRING:
        case TYPE_MULTIPOLYGON:
        case TYPE_MULTISOLID:
        case TYPE_GEOMETRYCOLLECTION: {
          const auto &gc = geom.as<GeometryCollection>();
          for (size_t i = 0; i < gc.numGeometries(); ++i) {
            Nef_polyhedron_3 temp = geometryToNef(gc.geometryN(i));
            if (i == 0) {
              result = temp;
            } else {
              result += temp;
            }
          }
          break;
        }
        default:
          throw std::runtime_error("Unsupported geometry type: " +
                                   geom.geometryType());
        }
      };

  process_geometry(g);
  return result;
}

// Helper function to convert Nef_polyhedron_3 to SFCGAL::Geometry
auto
nefToGeometry(const Nef_polyhedron_3 &nef) -> std::unique_ptr<Geometry>
{
  if (nef.is_empty()) {
    return std::make_unique<GeometryCollection>();
  }

  Polyhedron_3 poly;
  nef.convert_to_polyhedron(poly);

  if (poly.is_empty()) {
    return std::make_unique<GeometryCollection>();
  }

  return std::make_unique<PolyhedralSurface>(poly);
}

auto
minkowskiSum3D(const Geometry &gA, const Geometry &gB, NoValidityCheck /*unused*/)
    -> std::unique_ptr<Geometry>
{
  if (gA.isEmpty() || gB.isEmpty()) {
    return std::unique_ptr<Geometry>(new GeometryCollection());
  }

  Nef_polyhedron_3 nefA = geometryToNef(gA);
  Nef_polyhedron_3 nefB = geometryToNef(gB);

  if (nefA.is_empty() || nefB.is_empty()) {
    return std::unique_ptr<Geometry>(new GeometryCollection());
  }

  Nef_polyhedron_3 result = CGAL::minkowski_sum_3(nefA, nefB);

  if (result.is_empty()) {
    return std::unique_ptr<Geometry>(new GeometryCollection());
  }

  return nefToGeometry(result);
}

auto
minkowskiSum3D(const Geometry &gA, const Geometry &gB)
    -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(gA);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(gB);

  std::unique_ptr<Geometry> result(minkowskiSum3D(gA, gB, NoValidityCheck()));
  propagateValidityFlag(*result, true);
  return result;
}

} // namespace SFCGAL::algorithm
