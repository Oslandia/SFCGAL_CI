#include "SFCGAL/algorithm/lexicographicOrder.h"
#include "SFCGAL/algorithm/orientation.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Triangle.h"

#include <algorithm>
#include <limits>

namespace SFCGAL {
namespace transform {

struct LexicographicPointComparator {
  bool operator()(const Point &p1, const Point &p2) const
  {
    if (p1.x() != p2.x())
      return p1.x() < p2.x();
    if (p1.y() != p2.y())
      return p1.y() < p2.y();

    if (p1.is3D() && p2.is3D()) {
      if (p1.z() != p2.z())
        return p1.z() < p2.z();
    } else if (p1.is3D() && !p2.is3D()) {
      return false;
    } else if (!p1.is3D() && p2.is3D()) {
      return true;
    }

    if (p1.isMeasured() && p2.isMeasured()) {
      if (p1.m() != p2.m())
        return p1.m() < p2.m();
    } else if (p1.isMeasured() && !p2.isMeasured()) {
      return false;
    } else if (!p1.isMeasured() && p2.isMeasured()) {
      return true;
    }

    return false;
  }
};

std::unique_ptr<LineString>
lexicographicOrderRing(const LineString &ring)
{
  if (ring.numPoints() <= 1)
    return std::make_unique<LineString>(ring);

  std::vector<Point> points;
  for (size_t i = 0; i < ring.numPoints(); ++i) {
    points.push_back(ring.pointN(i));
  }

  auto minIt = std::min_element(points.begin(), points.end(),
                                LexicographicPointComparator());

  std::rotate(points.begin(), minIt, points.end());

  if (!(points.front() == points.back())) {
    points.push_back(points.front());
  }

  LineString tempRing;
  for (const auto &p : points) {
    tempRing.addPoint(p);
  }

  bool isCCW = SFCGAL::algorithm::isCounterClockWiseOriented(tempRing);

  if (isCCW) {
    std::reverse(points.begin(), points.end() - 1);
    points.back() = points.front();
  }

  auto newRing = std::make_unique<LineString>();
  for (const auto &p : points) {
    newRing->addPoint(p);
  }

  return newRing;
}

Point
getRepresentativePoint(const Geometry &g)
{
  switch (g.geometryTypeId()) {
  case TYPE_POINT:
    return g.as<Point>();
  case TYPE_LINESTRING:
    return g.as<LineString>().pointN(0);
  case TYPE_POLYGON:
    return g.as<Polygon>().exteriorRing().pointN(0);
  case TYPE_TRIANGLE:
    return g.as<Triangle>().vertex(0);
  case TYPE_SOLID:
    return getRepresentativePoint(g.as<Solid>().exteriorShell());
  case TYPE_POLYHEDRALSURFACE:
    if (g.as<PolyhedralSurface>().numPolygons() > 0) {
      return getRepresentativePoint(g.as<PolyhedralSurface>().polygonN(0));
    }
    return Point();
  default:
    if (g.isEmpty()) {
      return Point();
    }
    return getRepresentativePoint(g.geometryN(0));
  }
}

std::unique_ptr<Geometry>
lexicographicOrder(const Geometry &geometry)
{
  switch (geometry.geometryTypeId()) {
  case TYPE_POINT:
  case TYPE_LINESTRING:
    return std::unique_ptr<Geometry>(geometry.clone());
  case TYPE_MULTIPOINT: {
    const MultiPoint  &mp = geometry.as<MultiPoint>();
    std::vector<Point> points;
    for (size_t i = 0; i < mp.numGeometries(); ++i) {
      points.push_back(mp.geometryN(i).as<Point>());
    }
    std::sort(points.begin(), points.end(), LexicographicPointComparator());
    auto newMP = std::make_unique<MultiPoint>();
    for (const auto &p : points) {
      newMP->addGeometry(new Point(p));
    }
    return newMP;
  }

  case TYPE_MULTILINESTRING: {
    const MultiLineString &mls = geometry.as<MultiLineString>();
    std::vector<std::pair<Point, std::unique_ptr<LineString>>> sortedLines;
    for (size_t i = 0; i < mls.numGeometries(); ++i) {
      const LineString &line = mls.geometryN(i).as<LineString>();
      sortedLines.emplace_back(line.pointN(0),
                               std::make_unique<LineString>(line));
    }
    std::sort(sortedLines.begin(), sortedLines.end(),
              [](const auto &a, const auto &b) {
                return LexicographicPointComparator()(a.first, b.first);
              });

    auto newMLS = std::make_unique<MultiLineString>();
    for (auto &pair : sortedLines) {
      newMLS->addGeometry(pair.second.release());
    }
    return newMLS;
  }
  case TYPE_POLYGON: {
    const Polygon &poly = geometry.as<Polygon>();

    auto orderedExteriorRing = lexicographicOrderRing(poly.exteriorRing());

    std::vector<std::unique_ptr<LineString>> orderedInteriorRings;
    for (size_t i = 0; i < poly.numInteriorRings(); ++i) {
      const LineString &interiorRing = poly.interiorRingN(i);
      auto orderedRing = lexicographicOrderRing(interiorRing);
      orderedInteriorRings.push_back(std::move(orderedRing));
    }

    std::sort(orderedInteriorRings.begin(), orderedInteriorRings.end(),
              [](const std::unique_ptr<LineString> &a, const std::unique_ptr<LineString> &b) {
                Point repA = getRepresentativePoint(*a);
                Point repB = getRepresentativePoint(*b);
                return LexicographicPointComparator()(repA, repB);
              });

    auto newPoly = std::make_unique<Polygon>();

    newPoly->setExteriorRing(orderedExteriorRing.release());

    for (auto &ring : orderedInteriorRings) {
      newPoly->addInteriorRing(ring.release());
    }

    return newPoly;
  }
  case TYPE_TRIANGLE: {
    const Triangle &triangle = geometry.as<Triangle>();

    std::vector<Point> vertices = {triangle.vertex(0), triangle.vertex(1), triangle.vertex(2)};

    auto minIt = std::min_element(vertices.begin(), vertices.end(), LexicographicPointComparator());
    std::rotate(vertices.begin(), minIt, vertices.end());

    LineString tempRing;
    tempRing.addPoint(vertices[0]);
    tempRing.addPoint(vertices[1]);
    tempRing.addPoint(vertices[2]);
    tempRing.addPoint(vertices[0]);

    bool isCCW = SFCGAL::algorithm::isCounterClockWiseOriented(tempRing);
    if (isCCW) {
      std::reverse(vertices.begin() + 1, vertices.end());
    }

    auto newTriangle = std::make_unique<Triangle>(vertices[0], vertices[1], vertices[2]);

    return newTriangle;
  }
  case TYPE_MULTIPOLYGON: {
    const MultiPolygon &mp = geometry.as<MultiPolygon>();

    std::vector<std::pair<Point, std::unique_ptr<Polygon>>> sortedPolygons;

    for (size_t i = 0; i < mp.numGeometries(); ++i) {
      const Polygon &poly = mp.geometryN(i).as<Polygon>();

      Point repPoint = getRepresentativePoint(poly);

      auto orderedPolygon = lexicographicOrder(poly);

      sortedPolygons.emplace_back(repPoint, std::unique_ptr<Polygon>(dynamic_cast<Polygon*>(orderedPolygon.release())));
    }

    std::sort(sortedPolygons.begin(), sortedPolygons.end(),
              [](const auto &a, const auto &b) {
                return LexicographicPointComparator()(a.first, b.first);
              });

    auto newMP = std::make_unique<MultiPolygon>();
    for (auto &pair : sortedPolygons) {
      newMP->addGeometry(pair.second.release());
    }
    return newMP;
  }
  case TYPE_POLYHEDRALSURFACE:
  case TYPE_TRIANGULATEDSURFACE: {
    const PolyhedralSurface &ps = geometry.as<PolyhedralSurface>();

    std::vector<std::pair<Point, std::unique_ptr<Polygon>>> sortedPolygons;

    for (size_t i = 0; i < ps.numPolygons(); ++i) {
      const Polygon &poly = ps.polygonN(i);

      Point repPoint = getRepresentativePoint(poly);

      auto orderedPolygon = lexicographicOrder(poly);

      sortedPolygons.emplace_back(repPoint, std::unique_ptr<Polygon>(dynamic_cast<Polygon*>(orderedPolygon.release())));
    }

    std::sort(sortedPolygons.begin(), sortedPolygons.end(),
              [](const auto &a, const auto &b) {
                return LexicographicPointComparator()(a.first, b.first);
              });

    auto newPS = std::make_unique<PolyhedralSurface>();
    for (auto &pair : sortedPolygons) {
      newPS->addPolygon(pair.second.release());
    }
    return newPS;
  }
  case TYPE_SOLID: {
    const Solid &solid = geometry.as<Solid>();

    auto orderedExteriorShellGeom = lexicographicOrder(solid.exteriorShell());
    auto orderedExteriorShell = std::unique_ptr<PolyhedralSurface>(dynamic_cast<PolyhedralSurface*>(orderedExteriorShellGeom.release()));

    std::vector<std::unique_ptr<PolyhedralSurface>> orderedInteriorShells;
    for (size_t i = 0; i < solid.numInteriorShells(); ++i) {
      auto orderedInteriorShellGeom = lexicographicOrder(solid.interiorShellN(i));
      auto orderedInteriorShell = std::unique_ptr<PolyhedralSurface>(dynamic_cast<PolyhedralSurface*>(orderedInteriorShellGeom.release()));
      orderedInteriorShells.push_back(std::move(orderedInteriorShell));
    }

    // Créer un nouveau Solid avec les coques triées
    auto newSolid = std::make_unique<Solid>(orderedExteriorShell.release());
    for (auto &shell : orderedInteriorShells) {
      newSolid->addInteriorShell(shell.release());
    }
    return newSolid;
  }

  case TYPE_MULTISOLID: {
    const MultiSolid &ms = geometry.as<MultiSolid>();

    std::vector<std::pair<Point, std::unique_ptr<Solid>>> sortedSolids;

    for (size_t i = 0; i < ms.numGeometries(); ++i) {
      const Solid &solid = ms.geometryN(i).as<Solid>();

      Point repPoint = getRepresentativePoint(solid);

      auto orderedSolid = lexicographicOrder(solid);

      sortedSolids.emplace_back(repPoint, std::unique_ptr<Solid>(dynamic_cast<Solid*>(orderedSolid.release())));
    }

    // Trier les solides
    std::sort(sortedSolids.begin(), sortedSolids.end(),
              [](const auto &a, const auto &b) {
                return LexicographicPointComparator()(a.first, b.first);
              });

    // Créer un nouveau MultiSolid avec les solides triés
    auto newMS = std::make_unique<MultiSolid>();
    for (auto &pair : sortedSolids) {
      newMS->addGeometry(pair.second.release());
    }
    return newMS;
  }

  case TYPE_GEOMETRYCOLLECTION: {
    const GeometryCollection &gc = geometry.as<GeometryCollection>();

    std::vector<std::unique_ptr<Geometry>> orderedGeometries;
    for (size_t i = 0; i < gc.numGeometries(); ++i) {
      auto orderedGeom = lexicographicOrder(gc.geometryN(i));
      orderedGeometries.push_back(std::move(orderedGeom));
    }

    std::sort(orderedGeometries.begin(), orderedGeometries.end(),
              [](const std::unique_ptr<Geometry> &a, const std::unique_ptr<Geometry> &b) {
                if (a->geometryTypeId() != b->geometryTypeId()) {
                  return a->geometryTypeId() < b->geometryTypeId();
                } else {
                  Point repA = getRepresentativePoint(*a);
                  Point repB = getRepresentativePoint(*b);
                  return LexicographicPointComparator()(repA, repB);
                }
              });

    auto newGC = std::make_unique<GeometryCollection>();
    for (auto &geom : orderedGeometries) {
      newGC->addGeometry(geom.release());
    }
    return newGC;
  }

  default:
    return std::unique_ptr<Geometry>(geometry.clone());
  }
}

} // namespace transform
} // namespace SFCGAL
