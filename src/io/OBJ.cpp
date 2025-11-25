// Copyright (c) 2024-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/io/OBJ.h"
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
#include "SFCGAL/numeric.h"
#include <algorithm>
#include <fstream>
#include <functional>
#include <limits>
#include <sstream>
#include <vector>

namespace SFCGAL::io::OBJ {

/**
 * @brief Parsed OBJ data structure
 */
struct ObjData {
  std::vector<Point>               vertices; ///< Vertex positions
  std::vector<std::vector<size_t>> faces;    ///< Face vertex indices
  std::vector<std::vector<size_t>> lines;    ///< Line vertex indices
  std::vector<size_t>              points;   ///< Point vertex indices
};

/**
 * @brief Parse OBJ data from input stream
 *
 * @param[in] inOBJ Input stream containing OBJ data
 * @return Parsed OBJ data
 * @throws SFCGAL::Exception If parsing fails
 */
// NOLINTBEGIN(readability-function-cognitive-complexity)
auto
parseObjData(std::istream &inOBJ) -> ObjData
{
  ObjData obj_data;

  std::string line;
  while (std::getline(inOBJ, line)) {
    // Skip empty lines and comments
    if (line.empty() || line[0] == '#') {
      continue;
    }

    std::istringstream iss(line);
    std::string        type;
    iss >> type;

    if (type == "v") {
      // Vertex: v x y z [w]
      double x = 0.0;
      double y = 0.0;
      double z = 0.0;
      iss >> x >> y;
      if (!(iss >> z)) {
        z = 0.0; // Default to 2D
      }
      obj_data.vertices.emplace_back(x, y, z);
    } else if (type == "f") {
      // Face: f v1 v2 v3 ...
      std::vector<size_t> face;
      std::string         vertex_data;
      while (iss >> vertex_data) {
        // Parse vertex index (handle v/vt/vn format by taking only first part)
        size_t      slash_pos        = vertex_data.find('/');
        std::string vertex_index_str = vertex_data.substr(0, slash_pos);

        try {
          size_t vertex_index = std::stoul(vertex_index_str);
          if (vertex_index == 0) {
            BOOST_THROW_EXCEPTION(
                Exception("Invalid vertex index 0 in OBJ file"));
          }
          // Convert from 1-based to 0-based indexing
          face.push_back(vertex_index - 1);
        } catch (const std::invalid_argument &) {
          BOOST_THROW_EXCEPTION(
              Exception("Invalid face vertex index: " + vertex_index_str));
        } catch (const std::out_of_range &) {
          BOOST_THROW_EXCEPTION(
              Exception("Face vertex index out of range: " + vertex_index_str));
        }
      }
      if (face.size() < 3) {
        BOOST_THROW_EXCEPTION(Exception("Face must have at least 3 vertices"));
      }
      obj_data.faces.push_back(face);
    } else if (type == "l") {
      // Line: l v1 v2 ...
      std::vector<size_t> line_indices;
      size_t              vertex_index;
      while (iss >> vertex_index) {
        if (vertex_index == 0) {
          BOOST_THROW_EXCEPTION(
              Exception("Invalid vertex index 0 in OBJ file"));
        }
        line_indices.push_back(vertex_index - 1);
      }
      if (line_indices.size() < 2) {
        BOOST_THROW_EXCEPTION(Exception("Line must have at least 2 vertices"));
      }
      obj_data.lines.push_back(line_indices);
    } else if (type == "p") {
      // Point: p v1
      size_t vertex_index;
      if (iss >> vertex_index) {
        if (vertex_index == 0) {
          BOOST_THROW_EXCEPTION(
              Exception("Invalid vertex index 0 in OBJ file"));
        }
        obj_data.points.push_back(vertex_index - 1);
      }
    }
    // Ignore other OBJ elements (materials, textures, etc.)
  }

  return obj_data;
}
// NOLINTEND(readability-function-cognitive-complexity)

/**
 * @brief Create geometry from parsed OBJ data
 *
 * @param[in] obj_data Parsed OBJ data
 * @return Geometry created from OBJ data
 * @throws SFCGAL::Exception If geometry creation fails
 */
// NOLINTBEGIN(readability-function-cognitive-complexity)
auto
createGeometryFromObjData(const ObjData &obj_data) -> std::unique_ptr<Geometry>
{
  const auto &vertices = obj_data.vertices;
  const auto &faces    = obj_data.faces;
  const auto &lines    = obj_data.lines;
  const auto &points   = obj_data.points;

  // Create geometry based on what we found
  if (faces.empty() && lines.empty() && points.empty()) {
    BOOST_THROW_EXCEPTION(Exception("No geometry found in OBJ file"));
  }

  // If we have faces, create a TriangulatedSurface or PolyhedralSurface
  if (!faces.empty()) {
    // Check if all faces are triangles
    // NOTE: for performance we could use only PolyhedralSurface and avoid this
    // "if" statement however, we could prefer TIN in most case
    bool all_triangles = std::all_of(
        faces.begin(), faces.end(),
        [](const std::vector<size_t> &face) { return face.size() == 3; });

    if (all_triangles) {
      // Create TriangulatedSurface
      auto triangulated_surface = std::make_unique<TriangulatedSurface>();

      for (const auto &face : faces) {
        if (face[0] >= vertices.size() || face[1] >= vertices.size() ||
            face[2] >= vertices.size()) {
          BOOST_THROW_EXCEPTION(
              Exception("Face references invalid vertex index"));
        }

        auto triangle = std::make_unique<Triangle>(
            vertices[face[0]], vertices[face[1]], vertices[face[2]]);
        triangulated_surface->addPatch(std::move(triangle));
      }

      return triangulated_surface;
    }
    // Create PolyhedralSurface
    auto polyhedral_surface = std::make_unique<PolyhedralSurface>();

    for (const auto &face : faces) {
      auto ring = std::make_unique<LineString>();

      for (size_t vertex_idx : face) {
        if (vertex_idx >= vertices.size()) {
          BOOST_THROW_EXCEPTION(
              Exception("Face references invalid vertex index"));
        }
        ring->addPoint(vertices[vertex_idx]);
      }
      // Close the ring
      if (!face.empty()) {
        ring->addPoint(vertices[face[0]]);
      }

      auto polygon = std::make_unique<Polygon>(ring.release());
      polyhedral_surface->addPatch(std::move(polygon));
    }

    return polyhedral_surface;
  }

  // If we only have lines, create a MultiLineString
  if (!lines.empty()) {
    auto multilinestring = std::make_unique<MultiLineString>();

    for (const auto &line_indices : lines) {
      auto linestring = std::make_unique<LineString>();
      for (size_t vertex_idx : line_indices) {
        if (vertex_idx >= vertices.size()) {
          BOOST_THROW_EXCEPTION(
              Exception("Line references invalid vertex index"));
        }
        linestring->addPoint(vertices[vertex_idx]);
      }
      multilinestring->addGeometry(std::move(linestring));
    }

    return multilinestring;
  }

  // If we only have points, create a MultiPoint
  auto multipoint = std::make_unique<MultiPoint>();

  for (size_t vertex_idx : points) {
    if (vertex_idx >= vertices.size()) {
      BOOST_THROW_EXCEPTION(Exception("Point references invalid vertex index"));
    }
    multipoint->addGeometry(std::make_unique<Point>(vertices[vertex_idx]));
  }

  return multipoint;
}
// NOLINTEND(readability-function-cognitive-complexity)

// NOLINTBEGIN(readability-function-cognitive-complexity)
void
save(const Geometry &geom, std::ostream &out)
{
  std::vector<Point>               unique_points;
  std::vector<std::vector<size_t>> all_faces;
  std::vector<std::vector<size_t>> all_lines;
  std::vector<size_t>              all_points_indices;

  const auto epsilon_ft = SFCGAL::Kernel::FT(SFCGAL::EPSILON);

  // Simple and deterministic point deduplication
  auto find_or_add_point = [&](const Point &point) -> size_t {
    // Linear search with epsilon comparison
    for (size_t i = 0; i < unique_points.size(); ++i) {
      const auto &existing = unique_points[i];
      if (almostEqual(point.x(), existing.x(), epsilon_ft) &&
          almostEqual(point.y(), existing.y(), epsilon_ft) &&
          almostEqual(point.is3D() ? point.z() : SFCGAL::Kernel::FT(0),
                      existing.is3D() ? existing.z() : SFCGAL::Kernel::FT(0),
                      epsilon_ft)) {
        return i + 1; // 1-based indexing
      }
    }

    // Add new point
    unique_points.push_back(point);
    return unique_points.size(); // 1-based indexing
  };

  std::function<void(const Geometry &)> process_geometry =
      [&](const Geometry &geom) {
        switch (geom.geometryTypeId()) {
        case TYPE_POINT: {
          const auto &point = geom.as<Point>();
          size_t      idx   = find_or_add_point(point);
          all_points_indices.push_back(idx);
          break;
        }
        case TYPE_LINESTRING: {
          const auto         &linestring = geom.as<LineString>();
          std::vector<size_t> line;
          for (size_t i = 0; i < linestring.numPoints(); ++i) {
            size_t idx = find_or_add_point(linestring.pointN(i));
            line.push_back(idx);
          }
          all_lines.push_back(line);
          break;
        }
        case TYPE_TRIANGLE: {
          const auto         &triangle = geom.as<Triangle>();
          std::vector<size_t> face;
          for (int i = 0; i < 3; ++i) {
            size_t idx = find_or_add_point(triangle.vertex(i));
            face.push_back(idx);
          }
          all_faces.push_back(face);
          break;
        }
        case TYPE_POLYGON: {
          const auto         &polygon = geom.as<Polygon>();
          std::vector<size_t> face;
          for (size_t i = 0; i < polygon.exteriorRing().numPoints() - 1; ++i) {
            size_t idx = find_or_add_point(polygon.exteriorRing().pointN(i));
            face.push_back(idx);
          }
          all_faces.push_back(face);
          break;
        }
        case TYPE_TRIANGULATEDSURFACE: {
          const auto &tin = geom.as<TriangulatedSurface>();
          for (size_t i = 0; i < tin.numPatches(); ++i) {
            process_geometry(tin.patchN(i));
          }
          break;
        }
        case TYPE_POLYHEDRALSURFACE: {
          const auto &polyhedral_surface = geom.as<PolyhedralSurface>();
          for (size_t i = 0; i < polyhedral_surface.numPatches(); ++i) {
            process_geometry(polyhedral_surface.patchN(i));
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
          const auto &geometry_collection = geom.as<GeometryCollection>();
          for (size_t i = 0; i < geometry_collection.numGeometries(); ++i) {
            process_geometry(geometry_collection.geometryN(i));
          }
          break;
        }
        default:
          BOOST_THROW_EXCEPTION(InappropriateGeometryException(
              "Unsupported geometry type: " + geom.geometryType()));
        }
      };

  process_geometry(geom);

  // Use buffered output for better I/O performance with large meshes
  std::ostringstream buffer;
  buffer.precision(out.precision());

  // Write vertices
  for (const auto &point : unique_points) {
    buffer << "v " << point.x() << " " << point.y() << " "
           << (point.is3D() ? point.z() : SFCGAL::Kernel::FT(0)) << "\n";
  }

  // Write points
  for (size_t idx : all_points_indices) {
    buffer << "p " << idx << "\n";
  }

  // Write lines
  for (const auto &line : all_lines) {
    buffer << "l";
    for (size_t idx : line) {
      buffer << " " << idx;
    }
    buffer << "\n";
  }

  // Write faces
  for (const auto &face : all_faces) {
    buffer << "f";
    for (size_t idx : face) {
      buffer << " " << idx;
    }
    buffer << "\n";
  }

  // Write all buffered content at once
  out << buffer.str();
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

auto
load(std::istream &inOBJ) -> std::unique_ptr<Geometry>
{
  ObjData obj_data = parseObjData(inOBJ);
  return createGeometryFromObjData(obj_data);
}

auto
load(const std::string &obj) -> std::unique_ptr<Geometry>
{
  std::istringstream iss(obj);
  return load(iss);
}

auto
loadFromFile(const std::string &filename) -> std::unique_ptr<Geometry>
{
  std::ifstream inOBJ(filename);
  if (!inOBJ) {
    BOOST_THROW_EXCEPTION(
        Exception("Unable to open file " + filename + " for reading."));
  }
  return load(inOBJ);
}

} // namespace SFCGAL::io::OBJ
