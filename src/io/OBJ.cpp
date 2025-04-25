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
#include <fstream>
#include <functional>
#include <sstream>
#include <vector>

namespace SFCGAL::io::OBJ {

void
save(const Geometry &geom, std::ostream &out)
{
  std::vector<Point>               all_points;
  std::vector<std::vector<size_t>> all_faces;
  std::vector<std::vector<size_t>> all_lines;
  std::vector<size_t>              all_points_indices;

  std::function<void(const Geometry &)> process_geometry =
      [&](const Geometry &g) {
        switch (g.geometryTypeId()) {
        case TYPE_POINT: {
          const auto &p = g.as<Point>();
          all_points.push_back(p);
          all_points_indices.push_back(all_points.size());
          break;
        }
        case TYPE_LINESTRING: {
          const auto         &ls = g.as<LineString>();
          std::vector<size_t> line;
          for (size_t i = 0; i < ls.numPoints(); ++i) {
            all_points.push_back(ls.pointN(i));
            line.push_back(all_points.size());
          }
          all_lines.push_back(line);
          break;
        }
        case TYPE_TRIANGLE: {
          const auto         &tri = g.as<Triangle>();
          std::vector<size_t> face;
          for (int i = 0; i < 3; ++i) {
            all_points.push_back(tri.vertex(i));
            face.push_back(all_points.size());
          }
          all_faces.push_back(face);
          break;
        }
        case TYPE_POLYGON: {
          const auto         &poly = g.as<Polygon>();
          std::vector<size_t> face;
          for (size_t i = 0; i < poly.exteriorRing().numPoints() - 1; ++i) {
            all_points.push_back(poly.exteriorRing().pointN(i));
            face.push_back(all_points.size());
          }
          all_faces.push_back(face);
          break;
        }
        case TYPE_TRIANGULATEDSURFACE: {
          const auto &ts = g.as<TriangulatedSurface>();
          for (size_t i = 0; i < ts.numPatchs(); ++i) {
            process_geometry(ts.patchN(i));
          }
          break;
        }
        case TYPE_POLYHEDRALSURFACE: {
          const auto &ps = g.as<PolyhedralSurface>();
          for (size_t i = 0; i < ps.numPatchs(); ++i) {
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

  for (const auto &p : all_points) {
    out << "v " << p.x() << " " << p.y() << " " << (p.is3D() ? p.z() : 0.0)
        << "\n";
  }

  for (size_t idx : all_points_indices) {
    out << "p " << idx << "\n";
  }

  for (const auto &line : all_lines) {
    out << "l";
    for (size_t idx : line) {
      out << " " << idx;
    }
    out << "\n";
  }

  for (const auto &face : all_faces) {
    out << "f";
    for (size_t idx : face) {
      out << " " << idx;
    }
    out << "\n";
  }
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

} // namespace SFCGAL::io::OBJ
