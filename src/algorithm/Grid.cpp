// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/Grid.h"

#include "SFCGAL/Envelope.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/detail/EnvelopeVisitor.h"

#include <cmath>

namespace SFCGAL::algorithm {

namespace {

/**
 * @brief Compute the envelope of a geometry
 */
auto
computeEnvelope(const Geometry &geom) -> Envelope
{
  Envelope                envelope;
  detail::EnvelopeVisitor visitor(envelope);
  geom.accept(visitor);
  return envelope;
}

/**
 * @brief Create a rectangle polygon at the given origin
 */
auto
createRectangle(double x, double y, double width, double height)
    -> std::unique_ptr<Polygon>
{
  auto ring = std::make_unique<LineString>();
  ring->addPoint(std::make_unique<Point>(x, y));
  ring->addPoint(std::make_unique<Point>(x + width, y));
  ring->addPoint(std::make_unique<Point>(x + width, y + height));
  ring->addPoint(std::make_unique<Point>(x, y + height));
  ring->addPoint(std::make_unique<Point>(x, y)); // close ring
  return std::make_unique<Polygon>(ring.release());
}

/**
 * @brief Create a flat-top hexagon polygon centered at (cx, cy)
 *
 * Flat-top hexagon has horizontal edges at top and bottom.
 * Vertices are at angles: 0, 60, 120, 180, 240, 300 degrees
 */
auto
createHexagonFlatTop(double cx, double cy, double size)
    -> std::unique_ptr<Polygon>
{
  auto         ring  = std::make_unique<LineString>();
  const double sqrt3 = std::sqrt(3.0);

  // Vertices at angles 0, 60, 120, 180, 240, 300 degrees
  // For flat-top: vertex at angle 0 is at (cx + size, cy)
  ring->addPoint(std::make_unique<Point>(cx + size, cy));
  ring->addPoint(std::make_unique<Point>(cx + size / 2, cy + size * sqrt3 / 2));
  ring->addPoint(std::make_unique<Point>(cx - size / 2, cy + size * sqrt3 / 2));
  ring->addPoint(std::make_unique<Point>(cx - size, cy));
  ring->addPoint(std::make_unique<Point>(cx - size / 2, cy - size * sqrt3 / 2));
  ring->addPoint(std::make_unique<Point>(cx + size / 2, cy - size * sqrt3 / 2));
  ring->addPoint(std::make_unique<Point>(cx + size, cy)); // close ring

  return std::make_unique<Polygon>(ring.release());
}

/**
 * @brief Create a pointy-top hexagon polygon centered at (cx, cy)
 *
 * Pointy-top hexagon has vertices at top and bottom.
 * Vertices are at angles: 30, 90, 150, 210, 270, 330 degrees
 */
auto
createHexagonPointyTop(double cx, double cy, double size)
    -> std::unique_ptr<Polygon>
{
  auto         ring  = std::make_unique<LineString>();
  const double sqrt3 = std::sqrt(3.0);

  // Vertices at angles 30, 90, 150, 210, 270, 330 degrees
  ring->addPoint(std::make_unique<Point>(cx + size * sqrt3 / 2, cy + size / 2));
  ring->addPoint(std::make_unique<Point>(cx, cy + size));
  ring->addPoint(std::make_unique<Point>(cx - size * sqrt3 / 2, cy + size / 2));
  ring->addPoint(std::make_unique<Point>(cx - size * sqrt3 / 2, cy - size / 2));
  ring->addPoint(std::make_unique<Point>(cx, cy - size));
  ring->addPoint(std::make_unique<Point>(cx + size * sqrt3 / 2, cy - size / 2));
  ring->addPoint(
      std::make_unique<Point>(cx + size * sqrt3 / 2, cy + size / 2)); // close

  return std::make_unique<Polygon>(ring.release());
}

/**
 * @brief Create a diamond polygon centered at (cx, cy)
 *
 * Diamond is a rotated square with vertices on the cardinal directions.
 */
auto
createDiamond(double cx, double cy, double size) -> std::unique_ptr<Polygon>
{
  auto         ring     = std::make_unique<LineString>();
  const double halfDiag = size / 2.0;

  ring->addPoint(std::make_unique<Point>(cx, cy + halfDiag)); // top
  ring->addPoint(std::make_unique<Point>(cx + halfDiag, cy)); // right
  ring->addPoint(std::make_unique<Point>(cx, cy - halfDiag)); // bottom
  ring->addPoint(std::make_unique<Point>(cx - halfDiag, cy)); // left
  ring->addPoint(std::make_unique<Point>(cx, cy + halfDiag)); // close

  return std::make_unique<Polygon>(ring.release());
}

/**
 * @brief Create a Solid cube at the given origin
 */
auto
createCube(double x, double y, double z, double sizeX, double sizeY,
           double sizeZ) -> std::unique_ptr<Solid>
{
  // Create 8 corner points
  Point a(x, y, z);
  Point b(x + sizeX, y, z);
  Point c(x + sizeX, y + sizeY, z);
  Point d(x, y + sizeY, z);
  Point e(x, y, z + sizeZ);
  Point f(x + sizeX, y, z + sizeZ);
  Point g(x + sizeX, y + sizeY, z + sizeZ);
  Point h(x, y + sizeY, z + sizeZ);

  auto shell = std::make_unique<PolyhedralSurface>();

  // Bottom face: a, d, c, b (counter-clockwise from outside)
  {
    LineString ring;
    ring.addPoint(a);
    ring.addPoint(d);
    ring.addPoint(c);
    ring.addPoint(b);
    ring.addPoint(a);
    shell->addPatch(Polygon(ring));
  }

  // Top face: e, f, g, h
  {
    LineString ring;
    ring.addPoint(e);
    ring.addPoint(f);
    ring.addPoint(g);
    ring.addPoint(h);
    ring.addPoint(e);
    shell->addPatch(Polygon(ring));
  }

  // Front face: a, b, f, e
  {
    LineString ring;
    ring.addPoint(a);
    ring.addPoint(b);
    ring.addPoint(f);
    ring.addPoint(e);
    ring.addPoint(a);
    shell->addPatch(Polygon(ring));
  }

  // Back face: c, d, h, g
  {
    LineString ring;
    ring.addPoint(c);
    ring.addPoint(d);
    ring.addPoint(h);
    ring.addPoint(g);
    ring.addPoint(c);
    shell->addPatch(Polygon(ring));
  }

  // Right face: b, c, g, f
  {
    LineString ring;
    ring.addPoint(b);
    ring.addPoint(c);
    ring.addPoint(g);
    ring.addPoint(f);
    ring.addPoint(b);
    shell->addPatch(Polygon(ring));
  }

  // Left face: a, e, h, d
  {
    LineString ring;
    ring.addPoint(a);
    ring.addPoint(e);
    ring.addPoint(h);
    ring.addPoint(d);
    ring.addPoint(a);
    shell->addPatch(Polygon(ring));
  }

  return std::make_unique<Solid>(shell.release());
}

/**
 * @brief Check if a geometry intersects with the extent (for clipping)
 */
auto
intersectsExtent(const Geometry &cell, const Geometry &extent) -> bool
{
  // Use bounding box check first for efficiency
  Envelope cellEnv   = computeEnvelope(cell);
  Envelope extentEnv = computeEnvelope(extent);
  return Envelope::overlaps(cellEnv, extentEnv);
}

/**
 * @brief Clip a cell to the extent boundary
 */
auto
clipToExtent(const Geometry &cell, const Geometry &extent)
    -> std::unique_ptr<Geometry>
{
  return intersection(cell, extent);
}

} // anonymous namespace

//------------------------------------------------------------------------------
// 2D Grid Implementations
//------------------------------------------------------------------------------

auto
makeSquareGrid(const Geometry &extent, double cellSize, GridClipMode clipMode)
    -> Grid2D
{
  return makeRectangleGrid(extent, cellSize, cellSize, clipMode);
}

auto
makeRectangleGrid(const Geometry &extent, double cellSizeX, double cellSizeY,
                  GridClipMode clipMode) -> Grid2D
{
  if (cellSizeX <= 0 || cellSizeY <= 0) {
    BOOST_THROW_EXCEPTION(
        Exception("Cell size must be positive for rectangle grid"));
  }

  Envelope env = computeEnvelope(extent);
  if (env.isEmpty()) {
    return Grid2D();
  }

  const double minX = env.xMin();
  const double minY = env.yMin();
  const double maxX = env.xMax();
  const double maxY = env.yMax();

  // Calculate number of cells
  const auto nx = static_cast<int64_t>(std::ceil((maxX - minX) / cellSizeX));
  const auto ny = static_cast<int64_t>(std::ceil((maxY - minY) / cellSizeY));

  Grid2D grid;
  grid.reserve(static_cast<size_t>(nx * ny));

  for (int64_t j = 0; j < ny; ++j) {
    for (int64_t i = 0; i < nx; ++i) {
      double x = minX + static_cast<double>(i) * cellSizeX;
      double y = minY + static_cast<double>(j) * cellSizeY;

      auto cell = createRectangle(x, y, cellSizeX, cellSizeY);

      if (clipMode == GridClipMode::CLIP_TO_EXTENT) {
        auto clipped = clipToExtent(*cell, extent);
        if (clipped && !clipped->isEmpty()) {
          grid.emplace_back(std::move(clipped), i, j);
        }
      } else {
        // COVER_EXTENT: include all cells that overlap extent
        if (intersectsExtent(*cell, extent)) {
          grid.emplace_back(std::move(cell), i, j);
        }
      }
    }
  }

  return grid;
}

auto
makeHexagonGrid(const Geometry &extent, double cellSize, bool flatTop,
                GridClipMode clipMode) -> Grid2D
{
  if (cellSize <= 0) {
    BOOST_THROW_EXCEPTION(
        Exception("Cell size must be positive for hexagon grid"));
  }

  Envelope env = computeEnvelope(extent);
  if (env.isEmpty()) {
    return Grid2D();
  }

  const double minX  = env.xMin();
  const double minY  = env.yMin();
  const double maxX  = env.xMax();
  const double maxY  = env.yMax();
  const double sqrt3 = std::sqrt(3.0);

  Grid2D grid;

  if (flatTop) {
    // Flat-top hexagon grid
    // Horizontal spacing: 1.5 * cellSize
    // Vertical spacing: sqrt(3) * cellSize
    // Odd columns offset by sqrt(3)/2 * cellSize
    const double hSpacing = 1.5 * cellSize;
    const double vSpacing = sqrt3 * cellSize;

    const auto nx = static_cast<int64_t>(
        std::ceil((maxX - minX + cellSize) / hSpacing) + 1);
    const auto ny = static_cast<int64_t>(
        std::ceil((maxY - minY + vSpacing / 2) / vSpacing) + 1);

    for (int64_t i = 0; i < nx; ++i) {
      for (int64_t j = 0; j < ny; ++j) {
        double cx = minX + static_cast<double>(i) * hSpacing;
        double cy = minY + static_cast<double>(j) * vSpacing;

        // Offset odd columns
        if (i % 2 != 0) {
          cy += vSpacing / 2;
        }

        auto cell = createHexagonFlatTop(cx, cy, cellSize);

        if (clipMode == GridClipMode::CLIP_TO_EXTENT) {
          auto clipped = clipToExtent(*cell, extent);
          if (clipped && !clipped->isEmpty()) {
            grid.emplace_back(std::move(clipped), i, j);
          }
        } else {
          if (intersectsExtent(*cell, extent)) {
            grid.emplace_back(std::move(cell), i, j);
          }
        }
      }
    }
  } else {
    // Pointy-top hexagon grid
    // Horizontal spacing: sqrt(3) * cellSize
    // Vertical spacing: 1.5 * cellSize
    // Odd rows offset by sqrt(3)/2 * cellSize
    const double hSpacing = sqrt3 * cellSize;
    const double vSpacing = 1.5 * cellSize;

    const auto nx = static_cast<int64_t>(
        std::ceil((maxX - minX + hSpacing / 2) / hSpacing) + 1);
    const auto ny = static_cast<int64_t>(
        std::ceil((maxY - minY + cellSize) / vSpacing) + 1);

    for (int64_t j = 0; j < ny; ++j) {
      for (int64_t i = 0; i < nx; ++i) {
        double cx = minX + static_cast<double>(i) * hSpacing;
        double cy = minY + static_cast<double>(j) * vSpacing;

        // Offset odd rows
        if (j % 2 != 0) {
          cx += hSpacing / 2;
        }

        auto cell = createHexagonPointyTop(cx, cy, cellSize);

        if (clipMode == GridClipMode::CLIP_TO_EXTENT) {
          auto clipped = clipToExtent(*cell, extent);
          if (clipped && !clipped->isEmpty()) {
            grid.emplace_back(std::move(clipped), i, j);
          }
        } else {
          if (intersectsExtent(*cell, extent)) {
            grid.emplace_back(std::move(cell), i, j);
          }
        }
      }
    }
  }

  return grid;
}

auto
makeTriangleGrid(const Geometry &extent, double cellSize, GridClipMode clipMode)
    -> Grid2D
{
  if (cellSize <= 0) {
    BOOST_THROW_EXCEPTION(
        Exception("Cell size must be positive for triangle grid"));
  }

  Envelope env = computeEnvelope(extent);
  if (env.isEmpty()) {
    return Grid2D();
  }

  const double minX   = env.xMin();
  const double minY   = env.yMin();
  const double maxX   = env.xMax();
  const double maxY   = env.yMax();
  const double height = cellSize * std::sqrt(3.0) / 2.0;

  // Calculate number of cells
  const auto nx = static_cast<int64_t>(std::ceil((maxX - minX) / cellSize) + 1);
  const auto ny = static_cast<int64_t>(std::ceil((maxY - minY) / height) + 1);

  Grid2D grid;

  for (int64_t j = 0; j < ny; ++j) {
    double yBase = minY + static_cast<double>(j) * height;

    for (int64_t i = 0; i < nx * 2; ++i) {
      // Two triangles per square: up-pointing and down-pointing
      bool   upPointing = (i + j) % 2 == 0;
      double xBase      = minX + static_cast<double>(i / 2) * cellSize;

      // Offset for odd rows
      if (j % 2 != 0) {
        xBase -= cellSize / 2;
      }

      auto ring = std::make_unique<LineString>();

      if (upPointing) {
        // Up-pointing triangle
        ring->addPoint(std::make_unique<Point>(xBase, yBase));
        ring->addPoint(std::make_unique<Point>(xBase + cellSize, yBase));
        ring->addPoint(
            std::make_unique<Point>(xBase + cellSize / 2, yBase + height));
        ring->addPoint(std::make_unique<Point>(xBase, yBase));
      } else {
        // Down-pointing triangle
        ring->addPoint(std::make_unique<Point>(xBase + cellSize / 2, yBase));
        ring->addPoint(
            std::make_unique<Point>(xBase + cellSize, yBase + height));
        ring->addPoint(std::make_unique<Point>(xBase, yBase + height));
        ring->addPoint(std::make_unique<Point>(xBase + cellSize / 2, yBase));
      }

      auto cell = std::make_unique<Polygon>(ring.release());

      if (clipMode == GridClipMode::CLIP_TO_EXTENT) {
        auto clipped = clipToExtent(*cell, extent);
        if (clipped && !clipped->isEmpty()) {
          grid.emplace_back(std::move(clipped), i, j);
        }
      } else {
        if (intersectsExtent(*cell, extent)) {
          grid.emplace_back(std::move(cell), i, j);
        }
      }
    }
  }

  return grid;
}

auto
makeDiamondGrid(const Geometry &extent, double cellSize, GridClipMode clipMode)
    -> Grid2D
{
  if (cellSize <= 0) {
    BOOST_THROW_EXCEPTION(
        Exception("Cell size must be positive for diamond grid"));
  }

  Envelope env = computeEnvelope(extent);
  if (env.isEmpty()) {
    return Grid2D();
  }

  const double minX     = env.xMin();
  const double minY     = env.yMin();
  const double maxX     = env.xMax();
  const double maxY     = env.yMax();
  const double halfDiag = cellSize / 2.0;

  // Diamond grid spacing
  const double hSpacing = cellSize;
  const double vSpacing = cellSize;

  const auto nx = static_cast<int64_t>(std::ceil((maxX - minX) / hSpacing) + 1);
  const auto ny = static_cast<int64_t>(std::ceil((maxY - minY) / vSpacing) + 1);

  Grid2D grid;

  for (int64_t j = 0; j < ny; ++j) {
    for (int64_t i = 0; i < nx; ++i) {
      double cx = minX + static_cast<double>(i) * hSpacing + halfDiag;
      double cy = minY + static_cast<double>(j) * vSpacing + halfDiag;

      // Offset odd rows
      if (j % 2 != 0) {
        cx += halfDiag;
      }

      auto cell = createDiamond(cx, cy, cellSize);

      if (clipMode == GridClipMode::CLIP_TO_EXTENT) {
        auto clipped = clipToExtent(*cell, extent);
        if (clipped && !clipped->isEmpty()) {
          grid.emplace_back(std::move(clipped), i, j);
        }
      } else {
        if (intersectsExtent(*cell, extent)) {
          grid.emplace_back(std::move(cell), i, j);
        }
      }
    }
  }

  return grid;
}

//------------------------------------------------------------------------------
// 3D Grid Implementations
//------------------------------------------------------------------------------

auto
makeVoxelGrid(const Geometry &extent, double cellSizeX, double cellSizeY,
              double cellSizeZ, GridClipMode clipMode) -> Grid3D
{
  if (cellSizeX <= 0 || cellSizeY <= 0 || cellSizeZ <= 0) {
    BOOST_THROW_EXCEPTION(
        Exception("Cell sizes must be positive for voxel grid"));
  }

  Envelope env = computeEnvelope(extent);
  if (env.isEmpty() || !env.is3D()) {
    BOOST_THROW_EXCEPTION(Exception("Extent must be a valid 3D geometry"));
  }

  const double minX = env.xMin();
  const double minY = env.yMin();
  const double minZ = env.zMin();
  const double maxX = env.xMax();
  const double maxY = env.yMax();
  const double maxZ = env.zMax();

  const auto nx = static_cast<int64_t>(std::ceil((maxX - minX) / cellSizeX));
  const auto ny = static_cast<int64_t>(std::ceil((maxY - minY) / cellSizeY));
  const auto nz = static_cast<int64_t>(std::ceil((maxZ - minZ) / cellSizeZ));

  Grid3D grid;
  grid.reserve(static_cast<size_t>(nx * ny * nz));

  for (int64_t k = 0; k < nz; ++k) {
    for (int64_t j = 0; j < ny; ++j) {
      for (int64_t i = 0; i < nx; ++i) {
        double x = minX + static_cast<double>(i) * cellSizeX;
        double y = minY + static_cast<double>(j) * cellSizeY;
        double z = minZ + static_cast<double>(k) * cellSizeZ;

        auto cell = createCube(x, y, z, cellSizeX, cellSizeY, cellSizeZ);

        if (clipMode == GridClipMode::CLIP_TO_EXTENT) {
          // 3D intersection is more complex, for now just check overlap
          // Full 3D clipping would require intersection3D
          Envelope cellEnv = computeEnvelope(*cell);
          if (Envelope::overlaps(cellEnv, env)) {
            auto clipped = intersection3D(*cell, extent);
            if (clipped && !clipped->isEmpty()) {
              grid.emplace_back(std::move(clipped), i, j, k);
            }
          }
        } else {
          Envelope cellEnv = computeEnvelope(*cell);
          if (Envelope::overlaps(cellEnv, env)) {
            grid.emplace_back(std::move(cell), i, j, k);
          }
        }
      }
    }
  }

  return grid;
}

auto
makeTetrahedronGrid(const Geometry &extent, double cellSize,
                    GridClipMode clipMode) -> Grid3D
{
  if (cellSize <= 0) {
    BOOST_THROW_EXCEPTION(
        Exception("Cell size must be positive for tetrahedron grid"));
  }

  Envelope env = computeEnvelope(extent);
  if (env.isEmpty() || !env.is3D()) {
    BOOST_THROW_EXCEPTION(Exception("Extent must be a valid 3D geometry"));
  }

  const double minX = env.xMin();
  const double minY = env.yMin();
  const double minZ = env.zMin();
  const double maxX = env.xMax();
  const double maxY = env.yMax();
  const double maxZ = env.zMax();

  const auto nx = static_cast<int64_t>(std::ceil((maxX - minX) / cellSize));
  const auto ny = static_cast<int64_t>(std::ceil((maxY - minY) / cellSize));
  const auto nz = static_cast<int64_t>(std::ceil((maxZ - minZ) / cellSize));

  Grid3D grid;

  // 5 tetrahedra per cube (Freudenthal triangulation)
  for (int64_t k = 0; k < nz; ++k) {
    for (int64_t j = 0; j < ny; ++j) {
      for (int64_t i = 0; i < nx; ++i) {
        double x = minX + static_cast<double>(i) * cellSize;
        double y = minY + static_cast<double>(j) * cellSize;
        double z = minZ + static_cast<double>(k) * cellSize;

        // 8 vertices of the cube
        Point v0(x, y, z);
        Point v1(x + cellSize, y, z);
        Point v2(x + cellSize, y + cellSize, z);
        Point v3(x, y + cellSize, z);
        Point v4(x, y, z + cellSize);
        Point v5(x + cellSize, y, z + cellSize);
        Point v6(x + cellSize, y + cellSize, z + cellSize);
        Point v7(x, y + cellSize, z + cellSize);

        // 5 tetrahedra from the cube
        // Tetrahedra vertices (using right-hand rule for outward normals)
        std::vector<std::array<const Point *, 4>> tetras = {
            {&v0, &v1, &v3, &v4}, // tet 0
            {&v1, &v2, &v3, &v6}, // tet 1
            {&v1, &v4, &v5, &v6}, // tet 2
            {&v3, &v4, &v6, &v7}, // tet 3
            {&v1, &v3, &v4, &v6}  // tet 4 (center)
        };

        int64_t tetIdx = 0;
        for (const auto &tet : tetras) {
          // Create tetrahedron as Solid with 4 triangular faces
          auto shell = std::make_unique<PolyhedralSurface>();

          // Face 0: v0, v1, v2
          {
            LineString ring;
            ring.addPoint(*tet[0]);
            ring.addPoint(*tet[2]);
            ring.addPoint(*tet[1]);
            ring.addPoint(*tet[0]);
            shell->addPatch(Polygon(ring));
          }
          // Face 1: v0, v1, v3
          {
            LineString ring;
            ring.addPoint(*tet[0]);
            ring.addPoint(*tet[1]);
            ring.addPoint(*tet[3]);
            ring.addPoint(*tet[0]);
            shell->addPatch(Polygon(ring));
          }
          // Face 2: v0, v2, v3
          {
            LineString ring;
            ring.addPoint(*tet[0]);
            ring.addPoint(*tet[3]);
            ring.addPoint(*tet[2]);
            ring.addPoint(*tet[0]);
            shell->addPatch(Polygon(ring));
          }
          // Face 3: v1, v2, v3
          {
            LineString ring;
            ring.addPoint(*tet[1]);
            ring.addPoint(*tet[2]);
            ring.addPoint(*tet[3]);
            ring.addPoint(*tet[1]);
            shell->addPatch(Polygon(ring));
          }

          auto cell = std::make_unique<Solid>(shell.release());

          // Combine grid indices: i, j encodes cube position
          // k encodes z-level, we use a sub-index for tetrahedra
          int64_t gridI = i * 5 + tetIdx;

          Envelope cellEnv = computeEnvelope(*cell);
          if (Envelope::overlaps(cellEnv, env)) {
            if (clipMode == GridClipMode::CLIP_TO_EXTENT) {
              auto clipped = intersection3D(*cell, extent);
              if (clipped && !clipped->isEmpty()) {
                grid.emplace_back(std::move(clipped), gridI, j, k);
              }
            } else {
              grid.emplace_back(std::move(cell), gridI, j, k);
            }
          }
          ++tetIdx;
        }
      }
    }
  }

  return grid;
}

//------------------------------------------------------------------------------
// Utility Functions
//------------------------------------------------------------------------------

auto
gridToGeometryCollection(Grid2D &grid) -> std::unique_ptr<GeometryCollection>
{
  auto gc = std::make_unique<GeometryCollection>();
  for (auto &cell : grid) {
    if (cell.geometry) {
      gc->addGeometry(cell.geometry.release());
    }
  }
  return gc;
}

auto
gridToGeometryCollection(Grid3D &grid) -> std::unique_ptr<GeometryCollection>
{
  auto gc = std::make_unique<GeometryCollection>();
  for (auto &cell : grid) {
    if (cell.geometry) {
      gc->addGeometry(cell.geometry.release());
    }
  }
  return gc;
}

} // namespace SFCGAL::algorithm
