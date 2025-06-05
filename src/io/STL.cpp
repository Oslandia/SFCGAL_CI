// Copyright (c) 2025-2025, Oslandia.
// Copyright (c) 2025-2025, SFCGAL team.
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

// Utility function to normalize "-0" strings to "0" for consistent output
auto
normalizeZeroString(const std::string &str) -> std::string
{
  return (str == "-0") ? "0" : str;
}

// Convert value using kernel and normalize -0 to 0
template <typename T>
auto
normalizeKernelValue(const T &value) -> std::string
{
  std::ostringstream oss;
  oss << value;
  return normalizeZeroString(oss.str());
}

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

    // Get normalized string values to avoid -0 output
    std::string normalizedX = normalizeKernelValue(normal.x());
    std::string normalizedY = normalizeKernelValue(normal.y());
    std::string normalizedZ = normalizeKernelValue(normal.z());

    out << "  facet normal " << normalizedX << " " << normalizedY << " "
        << normalizedZ << "\n";
    out << "    outer loop\n";
    for (int i = 0; i < 3; ++i) {
      const auto &vertex  = triangle.vertex(i);
      std::string vertexX = normalizeKernelValue(vertex.x());
      std::string vertexY = normalizeKernelValue(vertex.y());
      std::string vertexZ =
          vertex.is3D() ? normalizeKernelValue(vertex.z()) : "0";
      out << "      vertex " << vertexX << " " << vertexY << " " << vertexZ
          << "\n";
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
