// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_GRID_H_
#define SFCGAL_ALGORITHM_GRID_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/config.h"

#include <cstdint>
#include <memory>
#include <vector>

namespace SFCGAL::algorithm {

/**
 * @brief Enumeration of 2D grid types
 */
enum class GridType {
  SQUARE,        ///< Regular square cells
  RECTANGLE,     ///< Rectangular cells (different width/height)
  HEXAGON_FLAT,  ///< Flat-top hexagonal cells
  HEXAGON_POINTY ///< Pointy-top hexagonal cells
};

/**
 * @brief Enumeration of grid clipping modes
 */
enum class GridClipMode {
  CLIP_TO_EXTENT, ///< Clip cells to extent boundary
  COVER_EXTENT    ///< Full cells covering extent (may extend beyond)
};

/**
 * @brief A 2D grid cell with geometry and grid indices
 */
struct SFCGAL_API GridCell2D {
  std::unique_ptr<Geometry> geometry; ///< The cell geometry (Polygon/Triangle)
  int64_t                   i;        ///< Column index
  int64_t                   j;        ///< Row index

  /**
   * @brief Default constructor
   */
  GridCell2D() : geometry(nullptr), i(0), j(0) {}

  /**
   * @brief Constructor with geometry and indices
   * @param geom The cell geometry
   * @param col Column index
   * @param row Row index
   */
  GridCell2D(std::unique_ptr<Geometry> geom, int64_t col, int64_t row)
      : geometry(std::move(geom)), i(col), j(row)
  {
  }

  /**
   * @brief Move constructor
   * @param other The other GridCell2D object to move from
   */
  GridCell2D(GridCell2D &&other) noexcept
      : geometry(std::move(other.geometry)), i(other.i), j(other.j)
  {
  }

  /**
   * @brief Move assignment operator
   * @param other The other GridCell2D object to move from
   * @return Reference to this GridCell2D object
   */
  auto
  operator=(GridCell2D &&other) noexcept -> GridCell2D &
  {
    if (this != &other) {
      geometry = std::move(other.geometry);
      i        = other.i;
      j        = other.j;
    }
    return *this;
  }

  /// @brief Copy constructor (deleted)
  GridCell2D(const GridCell2D &other) = delete;

  /// @brief Copy assignment operator (deleted)
  auto
  operator=(const GridCell2D &other) -> GridCell2D & = delete;
};

/**
 * @brief A 3D grid cell with geometry and grid indices
 */
struct SFCGAL_API GridCell3D {
  std::unique_ptr<Geometry> geometry; ///< The cell geometry (Solid)
  int64_t                   i;        ///< X index
  int64_t                   j;        ///< Y index
  int64_t                   k;        ///< Z index

  /**
   * @brief Default constructor
   */
  GridCell3D() : geometry(nullptr), i(0), j(0), k(0) {}

  /**
   * @brief Constructor with geometry and indices
   * @param geom The cell geometry
   * @param xi X index
   * @param yj Y index
   * @param zk Z index
   */
  GridCell3D(std::unique_ptr<Geometry> geom, int64_t xi, int64_t yj, int64_t zk)
      : geometry(std::move(geom)), i(xi), j(yj), k(zk)
  {
  }

  /**
   * @brief Move constructor
   * @param other The other GridCell3D object to move from
   */
  GridCell3D(GridCell3D &&other) noexcept
      : geometry(std::move(other.geometry)), i(other.i), j(other.j), k(other.k)
  {
  }

  /**
   * @brief Move assignment operator
   * @param other The other GridCell3D object to move from
   * @return Reference to this GridCell3D object
   */
  auto
  operator=(GridCell3D &&other) noexcept -> GridCell3D &
  {
    if (this != &other) {
      geometry = std::move(other.geometry);
      i        = other.i;
      j        = other.j;
      k        = other.k;
    }
    return *this;
  }

  /// @brief Copy constructor (deleted)
  GridCell3D(const GridCell3D &other) = delete;

  /// @brief Copy assignment operator (deleted)
  auto
  operator=(const GridCell3D &other) -> GridCell3D & = delete;
};

/// @brief Type alias for a 2D grid
using Grid2D = std::vector<GridCell2D>;

/// @brief Type alias for a 3D grid
using Grid3D = std::vector<GridCell3D>;

//------------------------------------------------------------------------------
// 2D Grid Generation Functions
//------------------------------------------------------------------------------

/**
 * @brief Generate a square grid covering an extent
 *
 * Creates a regular grid of square cells that covers the given extent.
 * Each cell is a Polygon with its grid indices (i, j).
 *
 * @param extent The extent polygon to cover
 * @param cellSize The side length of each square cell
 * @param clipMode How to handle cells at the extent boundary
 * @return A vector of GridCell2D containing the grid cells
 * @throws Exception if cellSize <= 0 or extent is invalid
 *
 * @code
 * auto extent = io::readWkt("POLYGON((0 0,100 0,100 100,0 100,0 0))");
 * auto grid = algorithm::makeSquareGrid(*extent, 10.0);
 * // Creates 100 cells (10x10)
 * @endcode
 */
SFCGAL_API auto
makeSquareGrid(const Geometry &extent, double cellSize,
               GridClipMode clipMode = GridClipMode::COVER_EXTENT) -> Grid2D;

/**
 * @brief Generate a rectangular grid covering an extent
 *
 * Creates a regular grid of rectangular cells with different width and height.
 *
 * @param extent The extent polygon to cover
 * @param cellSizeX The width of each cell
 * @param cellSizeY The height of each cell
 * @param clipMode How to handle cells at the extent boundary
 * @return A vector of GridCell2D containing the grid cells
 * @throws Exception if cellSizeX <= 0 or cellSizeY <= 0 or extent is invalid
 */
SFCGAL_API auto
makeRectangleGrid(const Geometry &extent, double cellSizeX, double cellSizeY,
                  GridClipMode clipMode = GridClipMode::COVER_EXTENT) -> Grid2D;

/**
 * @brief Generate a hexagonal grid covering an extent
 *
 * Creates a regular grid of hexagonal cells.
 *
 * @param extent The extent polygon to cover
 * @param cellSize The edge length of each hexagon
 * @param flatTop If true, hexagons have a flat top edge; if false, pointy top
 * @param clipMode How to handle cells at the extent boundary
 * @return A vector of GridCell2D containing the grid cells
 * @throws Exception if cellSize <= 0 or extent is invalid
 *
 * @code
 * auto extent = io::readWkt("POLYGON((0 0,100 0,100 80,0 80,0 0))");
 * auto grid = algorithm::makeHexagonGrid(*extent, 10.0, true);
 * @endcode
 */
SFCGAL_API auto
makeHexagonGrid(const Geometry &extent, double cellSize, bool flatTop = true,
                GridClipMode clipMode = GridClipMode::COVER_EXTENT) -> Grid2D;

/**
 * @brief Generate a triangular grid covering an extent
 *
 * Creates a regular grid of triangular cells.
 *
 * @param extent The extent polygon to cover
 * @param cellSize The edge length of each triangle
 * @param clipMode How to handle cells at the extent boundary
 * @return A vector of GridCell2D containing the grid cells
 * @throws Exception if cellSize <= 0 or extent is invalid
 */
SFCGAL_API auto
makeTriangleGrid(const Geometry &extent, double cellSize,
                 GridClipMode clipMode = GridClipMode::COVER_EXTENT) -> Grid2D;

/**
 * @brief Generate a diamond grid covering an extent
 *
 * Creates a regular grid of diamond (rotated square) cells.
 *
 * @param extent The extent polygon to cover
 * @param cellSize The diagonal length of each diamond
 * @param clipMode How to handle cells at the extent boundary
 * @return A vector of GridCell2D containing the grid cells
 * @throws Exception if cellSize <= 0 or extent is invalid
 */
SFCGAL_API auto
makeDiamondGrid(const Geometry &extent, double cellSize,
                GridClipMode clipMode = GridClipMode::COVER_EXTENT) -> Grid2D;

//------------------------------------------------------------------------------
// 3D Grid Generation Functions
//------------------------------------------------------------------------------

/**
 * @brief Generate a 3D voxel (cube) grid covering an extent
 *
 * Creates a regular 3D grid of cubic cells. Each cell is a Solid.
 *
 * @param extent The 3D extent (Solid or geometry with 3D envelope)
 * @param cellSizeX The X dimension of each voxel
 * @param cellSizeY The Y dimension of each voxel
 * @param cellSizeZ The Z dimension of each voxel
 * @param clipMode How to handle cells at the extent boundary
 * @return A vector of GridCell3D containing the grid cells
 * @throws Exception if any cellSize <= 0 or extent is invalid
 */
SFCGAL_API auto
makeVoxelGrid(const Geometry &extent, double cellSizeX, double cellSizeY,
              double       cellSizeZ,
              GridClipMode clipMode = GridClipMode::COVER_EXTENT) -> Grid3D;

/**
 * @brief Generate a 3D tetrahedral grid covering an extent
 *
 * Creates a 3D grid of tetrahedral cells. Uses 5 tetrahedra per cube
 * for space-filling tessellation.
 *
 * @param extent The 3D extent (Solid or geometry with 3D envelope)
 * @param cellSize The base cube size for tetrahedral subdivision
 * @param clipMode How to handle cells at the extent boundary
 * @return A vector of GridCell3D containing the grid cells
 * @throws Exception if cellSize <= 0 or extent is invalid
 */
SFCGAL_API auto
makeTetrahedronGrid(const Geometry &extent, double cellSize,
                    GridClipMode clipMode = GridClipMode::COVER_EXTENT)
    -> Grid3D;

//------------------------------------------------------------------------------
// Utility Functions
//------------------------------------------------------------------------------

/**
 * @brief Convert a 2D grid to a GeometryCollection
 *
 * Extracts the geometry from each cell and combines them into a
 * GeometryCollection. The cell indices are discarded.
 *
 * @param grid The 2D grid to convert
 * @return A GeometryCollection containing all cell geometries
 */
SFCGAL_API auto
gridToGeometryCollection(Grid2D &grid) -> std::unique_ptr<GeometryCollection>;

/**
 * @brief Convert a 3D grid to a GeometryCollection
 *
 * Extracts the geometry from each cell and combines them into a
 * GeometryCollection. The cell indices are discarded.
 *
 * @param grid The 3D grid to convert
 * @return A GeometryCollection containing all cell geometries
 */
SFCGAL_API auto
gridToGeometryCollection(Grid3D &grid) -> std::unique_ptr<GeometryCollection>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_GRID_H_
