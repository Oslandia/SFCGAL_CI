// Copyright (c) 2024-2024, Oslandia.
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

void
save(const Geometry &geom, std::ostream &out)
{
  std::vector<Triangle> all_triangles;

  std::function<void(const Geometry &)> process_geometry =
      [&](const Geometry &g) {
        switch (g.geometryTypeId()) {
        case TYPE_TRIANGLE: {
          all_triangles.push_back(g.as<Triangle>());
          break;
        }
        case TYPE_POLYGON: {
          const auto         &poly = g.as<Polygon>();
          TriangulatedSurface ts;
          triangulate::triangulatePolygon3D(poly, ts);
          for (size_t i = 0; i < ts.numTriangles(); ++i) {
            all_triangles.push_back(ts.triangleN(i));
          }
          break;
        }
        case TYPE_TRIANGULATEDSURFACE: {
          const auto &ts = g.as<TriangulatedSurface>();
          for (size_t i = 0; i < ts.numTriangles(); ++i) {
            all_triangles.push_back(ts.triangleN(i));
          }
          break;
        }
        case TYPE_POLYHEDRALSURFACE: {
          const auto &ps = g.as<PolyhedralSurface>();
          for (size_t i = 0; i < ps.numPolygons(); ++i) {
            process_geometry(ps.polygonN(i));
          }
          break;
        }
        case TYPE_SOLID: {
          const auto &solid = g.as<Solid>();
          process_geometry(solid.exteriorShell());
          break;
        }
        case TYPE_MULTIPOLYGON:
        case TYPE_MULTISOLID:
        case TYPE_GEOMETRYCOLLECTION: {
          const auto &gc = g.as<GeometryCollection>();
          for (size_t i = 0; i < gc.numGeometries(); ++i) {
            process_geometry(gc.geometryN(i));
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
      const auto &v = triangle.vertex(i);
      out << "      vertex " << v.x() << " " << v.y() << " "
          << (v.is3D() ? v.z() : 0.0) << "\n";
    }
    out << "    endloop\n";
    out << "  endfacet\n";
  }

  // Write STL footer
  out << "endsolid SFCGAL_export\n";
}

void
save(const Geometry &geom, const std::string &filename)
{
  std::ofstream out(filename);
  if (!out) {
    throw std::runtime_error("Unable to open file " + filename +
                             " for writing.");
  }
  save(geom, out);
}

std::string
saveToString(const Geometry &geom)
{
  std::ostringstream oss;
  save(geom, oss);
  return oss.str();
}

void
saveToBuffer(const Geometry &geom, char *buffer, size_t *size)
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
