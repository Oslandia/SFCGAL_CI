// Copyright (c) 2025-2025, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/io/STL.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include <fstream>
#include <functional>
#include <sstream>
#include <vector>

namespace SFCGAL::io::STL {

auto
save(const Geometry &geom, std::ostream &out) -> void
{
  std::vector<Triangle> all_triangles;

  std::function<void(const Geometry &)> process_geometry =
      [&](const Geometry &geom) {
        switch (geom.geometryTypeId()) {
        case TYPE_TRIANGLE: {
          all_triangles.push_back(geom.as<Triangle>());
          break;
        }
        case TYPE_POLYGON: {
          const auto         &poly = geom.as<Polygon>();
          TriangulatedSurface tin;
          triangulate::triangulatePolygon3D(poly, tin);
          for (size_t i = 0; i < tin.numTriangles(); ++i) {
            all_triangles.push_back(tin.triangleN(i));
          }
          break;
        }
        case TYPE_TRIANGULATEDSURFACE: {
          const auto &tin = geom.as<TriangulatedSurface>();
          for (size_t i = 0; i < tin.numTriangles(); ++i) {
            all_triangles.push_back(tin.triangleN(i));
          }
          break;
        }
        case TYPE_POLYHEDRALSURFACE: {
          const auto &phs = geom.as<PolyhedralSurface>();
          for (size_t i = 0; i < phs.numPolygons(); ++i) {
            process_geometry(phs.polygonN(i));
          }
          break;
        }
        case TYPE_SOLID: {
          const auto &solid = geom.as<Solid>();
          process_geometry(solid.exteriorShell());
          break;
        }
        case TYPE_MULTIPOLYGON:
        case TYPE_MULTISOLID:
        case TYPE_GEOMETRYCOLLECTION: {
          const auto &geomcoll = geom.as<GeometryCollection>();
          for (size_t i = 0; i < geomcoll.numGeometries(); ++i) {
            process_geometry(geomcoll.geometryN(i));
          }
          break;
        }
        default:
          // Ignore other geometry types as they can't be represented in STL
          break;
        }
      };

  process_geometry(geom);

  // Write STL header
  out << "solid SFCGAL_export\n";

  // Write triangles
  for (const auto &triangle : all_triangles) {
    CGAL::Vector_3<Kernel> normal = CGAL::normal(
        triangle.vertex(0).toPoint_3(), triangle.vertex(1).toPoint_3(),
        triangle.vertex(2).toPoint_3());
    out << "  facet normal " << normal.x() << " " << normal.y() << " "
        << normal.z() << "\n";
    out << "    outer loop\n";
    for (int i = 0; i < 3; ++i) {
      const auto &vertex = triangle.vertex(i);
      out << "      vertex " << vertex.x() << " " << vertex.y() << " "
          << (vertex.is3D() ? vertex.z() : 0.0) << "\n";
    }
    out << "    endloop\n";
    out << "  endfacet\n";
  }

  // Write STL footer
  out << "endsolid SFCGAL_export\n";
}

auto
save(const Geometry &geom, const std::string &filename) -> void
{
  std::ofstream out(filename);
  if (!out) {
    throw std::runtime_error("Unable to open file " + filename +
                             " for writing.");
  }
  save(geom, out);
}

auto
saveToString(const Geometry &geom) -> std::string
{
  std::ostringstream oss;
  save(geom, oss);
  return oss.str();
}

auto
saveToBuffer(const Geometry &geom, char *buffer, size_t *size) -> void
{
  std::string result = saveToString(geom);
  if ((buffer != nullptr) && *size >= result.size()) {
    std::copy(result.begin(), result.end(), buffer);
    *size = result.size();
  } else {
    *size = result.size();
  }
}

} // namespace SFCGAL::io::STL
