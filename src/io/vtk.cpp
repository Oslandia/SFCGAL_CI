// Copyright (c) 2024-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/io/vtk.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include <fstream>
#include <functional>
#include <sstream>
#include <vector>

namespace SFCGAL::io::VTK {

/**
 * Maximum recursion depth for processing nested geometry collections.
 * Prevents stack overflow attacks (CWE-674).
 */
static constexpr int MAX_RECURSION_DEPTH = 32;

// NOLINTBEGIN(readability-function-cognitive-complexity)
void
save(const Geometry &geom, std::ostream &out)
{
  std::vector<Point>               all_points;
  std::vector<std::vector<size_t>> all_cells;
  std::vector<int>                 cell_types;

  std::function<void(const Geometry &, int)> process_geometry =
      [&](const Geometry &geometry, int depth) {
        // Check recursion depth to prevent stack overflow (CWE-674)
        if (depth > MAX_RECURSION_DEPTH) {
          throw std::runtime_error(
              "VTK export: maximum recursion depth exceeded");
        }

        switch (geometry.geometryTypeId()) {
        case TYPE_POINT: {
          const auto &point = geometry.as<Point>();
          all_points.push_back(point);
          all_cells.push_back({all_points.size() - 1});
          cell_types.push_back(1); // VTK_VERTEX
          break;
        }
        case TYPE_LINESTRING: {
          const auto         &linestring = geometry.as<LineString>();
          std::vector<size_t> line;
          for (size_t i = 0; i < linestring.numPoints(); ++i) {
            all_points.push_back(linestring.pointN(i));
            line.push_back(all_points.size() - 1);
          }
          all_cells.push_back(line);
          cell_types.push_back(4); // VTK_POLY_LINE
          break;
        }
        case TYPE_TRIANGLE: {
          const auto         &triangle = geometry.as<Triangle>();
          std::vector<size_t> face;
          for (int i = 0; i < 3; ++i) {
            all_points.push_back(triangle.vertex(i));
            face.push_back(all_points.size() - 1);
          }
          all_cells.push_back(face);
          cell_types.push_back(5); // VTK_TRIANGLE
          break;
        }
        case TYPE_POLYGON: {
          const auto         &polygon = geometry.as<Polygon>();
          std::vector<size_t> face;
          for (size_t i = 0; i < polygon.exteriorRing().numPoints() - 1; ++i) {
            all_points.push_back(polygon.exteriorRing().pointN(i));
            face.push_back(all_points.size() - 1);
          }
          all_cells.push_back(face);
          cell_types.push_back(7); // VTK_POLYGON
          break;
        }
        case TYPE_TRIANGULATEDSURFACE: {
          const auto &triangulatedSurface = geometry.as<TriangulatedSurface>();
          for (size_t i = 0; i < triangulatedSurface.numPatches(); ++i) {
            process_geometry(triangulatedSurface.patchN(i), depth + 1);
          }
          break;
        }
        case TYPE_POLYHEDRALSURFACE: {
          const auto &polyhedralSurface = geometry.as<PolyhedralSurface>();
          for (size_t i = 0; i < polyhedralSurface.numPatches(); ++i) {
            process_geometry(polyhedralSurface.patchN(i), depth + 1);
          }
          break;
        }
        case TYPE_SOLID: {
          const auto &solid = geometry.as<Solid>();
          if (!solid.isEmpty()) {
            process_geometry(solid.exteriorShell(), depth + 1);
          }
          break;
        }
        case TYPE_MULTIPOINT:
        case TYPE_MULTILINESTRING:
        case TYPE_MULTIPOLYGON:
        case TYPE_MULTISOLID:
        case TYPE_GEOMETRYCOLLECTION: {
          const auto &collection = geometry.as<GeometryCollection>();
          for (size_t i = 0; i < collection.numGeometries(); ++i) {
            process_geometry(collection.geometryN(i), depth + 1);
          }
          break;
        }
        default:
          BOOST_THROW_EXCEPTION(InappropriateGeometryException(
              "Unsupported geometry type: " + geometry.geometryType()));
        }
      };

  process_geometry(geom, 0);

  // Write VTK header
  out << "# vtk DataFile Version 2.0\n";
  out << "SFCGAL Geometry\n";
  out << "ASCII\n";
  out << "DATASET UNSTRUCTURED_GRID\n";

  // Write points
  out << "POINTS " << all_points.size() << " float\n";
  for (const auto &point : all_points) {
    out << point.x() << " " << point.y() << " "
        << (point.is3D() ? point.z() : 0.0) << "\n";
  }

  // Write cells
  size_t total_cell_size = 0;
  for (const auto &cell : all_cells) {
    total_cell_size += cell.size() + 1; // +1 for the size prefix
  }
  out << "CELLS " << all_cells.size() << " " << total_cell_size << "\n";
  for (const auto &cell : all_cells) {
    out << cell.size();
    for (size_t idx : cell) {
      out << " " << idx;
    }
    out << "\n";
  }

  // Write cell types
  out << "CELL_TYPES " << cell_types.size() << "\n";
  for (int type : cell_types) {
    out << type << "\n";
  }
}
// NOLINTEND(readability-function-cognitive-complexity)

void
save(const Geometry &geom, const std::string &filename)
{
  std::ofstream out(filename);
  if (!out) {
    BOOST_THROW_EXCEPTION(
        Exception("Unable to open file " + filename + " for writing."));
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

} // namespace SFCGAL::io::VTK
