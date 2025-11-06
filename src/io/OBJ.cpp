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

void
save(const Geometry &geom, std::ostream &out)
{
  std::vector<Point>               unique_points;
  std::vector<std::vector<size_t>> all_faces;
  std::vector<std::vector<size_t>> all_lines;
  std::vector<size_t>              all_points_indices;

  const auto epsilon_ft = SFCGAL::Kernel::FT(SFCGAL::EPSILON);

  // Use spatial hashing for O(1) average lookup instead of O(n) linear search
  std::map<std::tuple<long, long, long>, size_t> point_hash_map;

  // Hash function that discretizes coordinates for epsilon-based comparison
  auto hash_point = [&](const Point &p) -> std::tuple<long, long, long> {
    constexpr double scale = 1e9; // Scale factor for hashing
    auto x_hash = static_cast<long>(CGAL::to_double(p.x()) * scale);
    auto y_hash = static_cast<long>(CGAL::to_double(p.y()) * scale);
    auto z_hash = static_cast<long>(CGAL::to_double(p.is3D() ? p.z() : SFCGAL::Kernel::FT(0)) * scale);
    return std::make_tuple(x_hash, y_hash, z_hash);
  };

  auto find_or_add_point = [&](const Point &p) -> size_t {
    auto hash = hash_point(p);

    // Check if we have points in nearby hash buckets (for epsilon tolerance)
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          auto nearby_hash = std::make_tuple(
            std::get<0>(hash) + dx,
            std::get<1>(hash) + dy,
            std::get<2>(hash) + dz
          );

          auto it = point_hash_map.find(nearby_hash);
          if (it != point_hash_map.end()) {
            const auto &existing = unique_points[it->second - 1]; // Convert from 1-based

            if (almostEqual(p.x(), existing.x(), epsilon_ft) &&
                almostEqual(p.y(), existing.y(), epsilon_ft) &&
                almostEqual(p.is3D() ? p.z() : SFCGAL::Kernel::FT(0),
                           existing.is3D() ? existing.z() : SFCGAL::Kernel::FT(0), epsilon_ft)) {
              return it->second; // Already 1-based
            }
          }
        }
      }
    }

    // Add new point
    unique_points.push_back(p);
    size_t index = unique_points.size(); // 1-based indexing
    point_hash_map[hash] = index;
    return index;
  };

  std::function<void(const Geometry &)> process_geometry =
      [&](const Geometry &g) {
        switch (g.geometryTypeId()) {
        case TYPE_POINT: {
          const auto &p = g.as<Point>();
          size_t idx = find_or_add_point(p);
          all_points_indices.push_back(idx);
          break;
        }
        case TYPE_LINESTRING: {
          const auto         &ls = g.as<LineString>();
          std::vector<size_t> line;
          for (size_t i = 0; i < ls.numPoints(); ++i) {
            size_t idx = find_or_add_point(ls.pointN(i));
            line.push_back(idx);
          }
          all_lines.push_back(line);
          break;
        }
        case TYPE_TRIANGLE: {
          const auto         &tri = g.as<Triangle>();
          std::vector<size_t> face;
          for (int i = 0; i < 3; ++i) {
            size_t idx = find_or_add_point(tri.vertex(i));
            face.push_back(idx);
          }
          all_faces.push_back(face);
          break;
        }
        case TYPE_POLYGON: {
          const auto         &poly = g.as<Polygon>();
          std::vector<size_t> face;
          for (size_t i = 0; i < poly.exteriorRing().numPoints() - 1; ++i) {
            size_t idx = find_or_add_point(poly.exteriorRing().pointN(i));
            face.push_back(idx);
          }
          all_faces.push_back(face);
          break;
        }
        case TYPE_TRIANGULATEDSURFACE: {
          const auto &ts = g.as<TriangulatedSurface>();
          for (size_t i = 0; i < ts.numPatches(); ++i) {
            process_geometry(ts.patchN(i));
          }
          break;
        }
        case TYPE_POLYHEDRALSURFACE: {
          const auto &ps = g.as<PolyhedralSurface>();
          for (size_t i = 0; i < ps.numPatches(); ++i) {
            process_geometry(ps.patchN(i));
          }
          break;
        }
        case TYPE_SOLID: {
          const auto &solid = g.as<Solid>();
          process_geometry(solid.exteriorShell());
          break;
        }
        case TYPE_MULTIPOINT:
        case TYPE_MULTILINESTRING:
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
          BOOST_THROW_EXCEPTION(InappropriateGeometryException(
              "Unsupported geometry type: " + g.geometryType()));
        }
      };

  process_geometry(geom);

  // Use buffered output for better I/O performance with large meshes
  std::ostringstream buffer;
  buffer.precision(out.precision());

  // Write vertices
  for (const auto &p : unique_points) {
    buffer << "v " << p.x() << " " << p.y() << " " << (p.is3D() ? p.z() : 0.0)
           << "\n";
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
load(std::istream &in) -> std::unique_ptr<Geometry>
{
  std::vector<Point>               vertices;
  std::vector<std::vector<size_t>> faces;
  std::vector<std::vector<size_t>> lines;
  std::vector<size_t>              points;

  std::string line;
  while (std::getline(in, line)) {
    // Skip empty lines and comments
    if (line.empty() || line[0] == '#') {
      continue;
    }

    std::istringstream iss(line);
    std::string        type;
    iss >> type;

    if (type == "v") {
      // Vertex: v x y z [w]
      double x, y, z = 0.0;
      iss >> x >> y;
      if (!(iss >> z)) {
        z = 0.0; // Default to 2D
      }
      vertices.emplace_back(x, y, z);
    } else if (type == "f") {
      // Face: f v1 v2 v3 ...
      std::vector<size_t> face;
      std::string         vertex_data;
      while (iss >> vertex_data) {
        // Parse vertex index (handle v/vt/vn format by taking only first part)
        size_t      slash_pos        = vertex_data.find('/');
        std::string vertex_index_str = vertex_data.substr(0, slash_pos);

        size_t vertex_index = std::stoul(vertex_index_str);
        if (vertex_index == 0) {
          BOOST_THROW_EXCEPTION(
              Exception("Invalid vertex index 0 in OBJ file"));
        }
        // Convert from 1-based to 0-based indexing
        face.push_back(vertex_index - 1);
      }
      if (face.size() < 3) {
        BOOST_THROW_EXCEPTION(Exception("Face must have at least 3 vertices"));
      }
      faces.push_back(face);
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
      lines.push_back(line_indices);
    } else if (type == "p") {
      // Point: p v1
      size_t vertex_index;
      if (iss >> vertex_index) {
        if (vertex_index == 0) {
          BOOST_THROW_EXCEPTION(
              Exception("Invalid vertex index 0 in OBJ file"));
        }
        points.push_back(vertex_index - 1);
      }
    }
    // Ignore other OBJ elements (materials, textures, etc.)
  }

  // Create geometry based on what we found
  if (faces.empty() && lines.empty() && points.empty()) {
    BOOST_THROW_EXCEPTION(Exception("No geometry found in OBJ file"));
  }

  // If we have faces, create a TriangulatedSurface or PolyhedralSurface
  if (!faces.empty()) {
    // Check if all faces are triangles
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

      return std::move(triangulated_surface);
    } else {
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

      return std::move(polyhedral_surface);
    }
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

    return std::move(multilinestring);
  }

  // If we only have points, create a MultiPoint
  if (!points.empty()) {
    auto multipoint = std::make_unique<MultiPoint>();

    for (size_t vertex_idx : points) {
      if (vertex_idx >= vertices.size()) {
        BOOST_THROW_EXCEPTION(
            Exception("Point references invalid vertex index"));
      }
      multipoint->addGeometry(std::make_unique<Point>(vertices[vertex_idx]));
    }

    return std::move(multipoint);
  }

  BOOST_THROW_EXCEPTION(Exception("No valid geometry found in OBJ file"));
}

auto
load(const std::string &filename) -> std::unique_ptr<Geometry>
{
  std::ifstream in(filename);
  if (!in) {
    BOOST_THROW_EXCEPTION(
        Exception("Unable to open file " + filename + " for reading."));
  }
  return load(in);
}

auto
loadFromString(const std::string &obj) -> std::unique_ptr<Geometry>
{
  std::istringstream iss(obj);
  return load(iss);
}

} // namespace SFCGAL::io::OBJ
