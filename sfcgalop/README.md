#SFCGALOP - SFCGAL Geometry Operations CLI

A command-line interface for performing geometric operations using the SFCGAL library.

## Features

- **50+ Geometric Operations**: Comprehensive support for 2D and 3D operations
- **Terminal UI**: Unicode box-drawing tables with color support
- **Multiple Input/Output Formats**: WKT, WKB support
- **Robust Error Handling**: Detailed validation messages using SFCGAL's validation engine
- **Cross-Platform**: Works on Linux, macOS, FreeBSD, and Windows

## Building

### Requirements

- C++17 compatible compiler
- CMake 3.10+
- SFCGAL library (>= 1.5.0)
- CGAL (>= 5.6)
- Boost (>= 1.69)
- GMP and MPFR libraries

### Build Instructions

```bash
mkdir build
cd build
cmake ..
make -j8
```

## Usage

### Basic Syntax

```bash
sfcgalop [options] -a <WKT/WKB> [-b <WKT/WKB>] [operation] [params]
```

**Note**: If no operation is specified, the input geometry will be displayed in the requested format (useful for format conversion).

### Options

- `-a, --geom-a <WKT/WKB>` : First geometry (required, use "stdin" to read from stdin)
- `-b, --geom-b <WKT/WKB>` : Second geometry (for binary operations, use "stdin" to read from stdin)
- `--validate` : Validate input geometries before operation
- `-f, --format <fmt>` : Output format (wkt, wkb, txt/ewkt, obj, stl, vtk)
- `--precision <n>` : Output precision (decimal places)
- `--list` : List all available operations
- `--help` : Show help message
- `--help-op <op>` : Detailed help for specific operation
- `-t, --time` : Show execution time
- `-v, --verbose` : Verbose output

### Examples

#### Calculate area of a polygon

```bash
sfcgalop -a "POLYGON((0 0, 0 10, 10 10, 10 0, 0 0))" area
#Output : 100
```

#### Calculate distance between two points

```bash
sfcgalop -a "POINT(0 0)" -b "POINT(3 4)" distance
#Output : 5
```

#### Compute intersection with validation

```bash
sfcgalop --validate \
  -a "POLYGON((0 0, 4 0, 4 4, 0 4, 0 0))" \
  -b "POLYGON((2 2, 6 2, 6 6, 2 6, 2 2))" \
  intersection
#Output : POLYGON((2 2, 2 4, 4 4, 4 2, 2 2))
```

#### Read from stdin

```bash
echo "POLYGON((0 0, 0 1, 1 1, 1 0, 0 0))" | \
  sfcgalop -a stdin convexhull
```

#### Convert geometry format

```bash
sfcgalop -a "POLYGON((0 0, 0 10, 10 10, 10 0, 0 0))" -f wkb
#Output : Binary WKB representation
```

#### Convert geometry to OBJ format

```bash
sfcgalop -a "TRIANGLE((0 0 0, 1 0 0, 0 1 0, 0 0 0))" -f obj
```

#### Convert geometry to STL format

```bash
sfcgalop -a "SOLID((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)), ((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0)), ((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0)), ((0 0 1,0 1 1,0 1 0,0 0 0,0 0 1)), ((1 0 1,1 1 1,0 1 1,0 0 1,1 0 1)), ((1 0 0,1 0 1,0 0 1,0 0 0,1 0 0))))" -f stl
```

#### Convert geometry to VTK format

```bash
sfcgalop -a "TIN Z (((0 0 0, 0 0 1, 0 1 0, 0 0 0)), ((0 0 0, 0 1 0, 1 0 0, 0 0 0)))" -f vtk
```

#### Display geometry from file

```bash
sfcgalop -a geometry.wkt
#Output : Geometry displayed in WKT format
```

#### 3D operations

```bash
sfcgalop -a "SOLID((((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),
  ((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0)),
  ((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),
  ((1 1 1, 0 1 1, 0 0 1, 1 0 1, 1 1 1)),
  ((1 1 1, 1 0 1, 1 0 0, 1 1 0, 1 1 1)),
  ((1 1 1, 1 1 0, 0 1 0, 0 1 1, 1 1 1))))" volume
#Output : 1
```

## Available Operations

### Metrics
- `area` - Calculate the 2D area of a geometry
- `area3d` - Calculate the 3D surface area of a geometry
- `volume` - Calculate the 3D volume of a solid geometry
- `length` - Calculate the 2D length of linear geometries
- `length3d` - Calculate the 3D length of linear geometries
- `distance` - Calculate the 2D minimum distance between two geometries
- `distance3d` - Calculate the 3D minimum distance between two geometries

### Predicates
- `intersects` - Test if two geometries intersect in 2D
- `intersects3d` - Test if two geometries intersect in 3D
- `covers` - Test if geometry A completely covers geometry B
- `is_valid` - Test if geometry is topologically valid
- `is_simple` - Test if geometry has no self-intersections
- `is_closed` - Test if linear geometry forms a closed ring
- `is_3d` - Test if geometry has Z coordinates
- `is_measured` - Test if geometry has measure (M) coordinates
- `is_empty` - Test if geometry contains no points

### Set Operations
- `intersection` - Compute the geometric intersection of two geometries
- `intersection3d` - Compute the 3D geometric intersection of two geometries
- `difference` - Compute geometry A minus geometry B
- `difference3d` - Compute 3D geometry A minus geometry B
- `union` - Compute the geometric union of two geometries
- `union3d` - Compute the 3D geometric union of two geometries

### Construction
- `boundary` - Compute the topological boundary of a geometry
- `envelope` - Compute the minimum bounding rectangle
- `convexhull` - Compute the 2D convex hull of a geometry
- `convexhull3d` - Compute the 3D convex hull of a geometry
- `centroid` - Compute the geometric centroid of a geometry
- `straightskeleton` - Compute the straight skeleton of a polygon
- `extrude` - Extrude a 2D geometry to create a 3D solid
- `tesselate` - Tesselate a geometry into triangular faces
- `triangulate` - Triangulate a geometry (alias for tesselate)
- `offset` - Create an offset polygon at specified distance
- `buffer3d` - Create a 3D buffer around points and lines
- `minkowskisum` - Compute Minkowski sum of two geometries
- `minkowskisum3d` - Compute 3D Minkowski sum of two geometries
- `alphashapes` - Compute alpha shapes from point cloud
- `alphawrapping3d` - Create 3D alpha wrapping surface from points
- `linesubstring` - Extract substring from linestring by fraction

### Transformations
- `translate` - Translate geometry by specified offset (params: dx dy [dz])
- `rotate` - Rotate geometry around Z-axis by angle (params: angle [cx cy])
- `scale` - Scale geometry by specified factors (params: sx [sy [sz]])
- `force2d` - Remove Z coordinates to create 2D geometry
- `force3d` - Add Z coordinates to create 3D geometry (params: [z_value])
- `forcemeasured` - Add measure coordinates to geometry (params: [m_value])
- `simplify` - Simplify geometry by removing vertices within tolerance (params: tolerance)
- `polygonRepair` - Repairs invalid polygons using repair's method (params: method)

### Collections
- `collect` - Combine two geometries into a collection
- `collection_extract` - Extract polygons from geometry collection
- `collection_homogenize` - Convert collection to appropriate multi-type
- `collection_to_multi` - Convert collection to multi-geometry type

### Analysis
- `orientation` - Determine polygon ring orientation (clockwise/counter-clockwise)
- `visibility` - Compute visibility polygon from a point in polygon (params: px py)
- `partition` - Partition polygon into simpler pieces (params: y|x|g|o|2|a)
- `normal` - Compute surface normal vector for polygon/triangle

## Architecture

### Core Components

- **main.cpp**: Entry point and CLI argument parsing
- **operations/operations.cpp**: All geometric operations implementation
- **io.cpp**: Geometry I/O handling (WKT, WKB)
- **error_handler.cpp**: Validation and error reporting
- **text_ui.cpp**: Terminal UI with tables and formatting

### Design Principles

1. **Robust Error Handling**: Detailed validation messages
2. **Extensible Architecture**: Easy to add new operations
3. **Cross-Platform Support**: Platform-specific handling isolated

## Error Handling

The tool provides detailed error messages with context:

```bash
#Invalid geometry
sfcgalop --validate -a "POLYGON((0 0, 0 10, 10 10, 10 0, 1 0))" area
#Output:
# ⚠ Geometry A validation issues:
#Validation Error : ring 0 is not closed
#Details:
#Geometry type : Polygon
#WKT : POLYGON((0 0, 0 10, 10 10, 10 0, 1 0))
```

## Naming Conventions and Compatibility

SFCGAL supports multiple naming conventions for better compatibility and flexibility:

### Underscore Convention Aliases
SFCGAL supports the `xx_yyy_zzz` naming convention as aliases for all operations:

- `line_substring` → `linesubstring`
- `alpha_shapes` → `alphashapes`
- `alpha_wrapping3d` → `alphawrapping3d`
- `alpha_wrapping_3d` → `alphawrapping3d`
- `alphawrapping_3d` → `alphawrapping3d`
- `buffer_3d` → `buffer3d`
- `minkowski_sum` → `minkowskisum`
- `minkowski_sum3d` → `minkowskisum3d`
- `minkowski_sum_3d` → `minkowskisum3d`
- `convex_hull` → `convexhull`
- `convex_hull3d` → `convexhull3d`
- `convex_hull_3d` → `convexhull3d`
- `straight_skeleton` → `straightskeleton`
- `force_2d` → `force2d`
- `force_3d` → `force3d`
- `force_measured` → `forcemeasured`
- `force_lhr` → `forceLHR`
- `force_rhr` → `forceRHR`
- `polygon_repair` → `polygonrepair`
- `distance_3d` → `distance3d`
- `length_3d` → `length3d`
- `area_3d` → `area3d`
- And many more

### GEOS Compatibility
SFCGAL also supports GEOS-style operation names as aliases:

- `buffer` → `offset` (GEOS name for offset operation)
- `symdifference` → `difference` (GEOS symmetric difference - note: SFCGAL uses regular difference)
- `geomunion` → `union` (Alternative to union)
- `geomintersection` → `intersection` (Alternative to intersection)
- `geomdifference` → `difference` (Alternative to difference)
- `isempty` → `is_empty` (Boolean check)
- `issimple` → `is_simple` (Boolean check)
- `isvalid` → `is_valid` (Boolean check)
- `is3d` → `is_3d` (Boolean check)
- `ismeasured` → `is_measured` (Boolean check)
- `isclosed` → `is_closed` (Boolean check)
- And many more common aliases

### Case Insensitivity
All operation names are case-insensitive, so `Area`, `AREA`, and `area` all work the same.

Example usage:
```bash
# Using underscore convention alias
sfcgalop -a "POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))" line_substring "start=0.25,end=0.75"
# Using GEOS-style alias
sfcgalop -a "POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))" isvalid
# Output: true
# Using original SFCGAL name
sfcgalop -a "POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))" is_valid
# Output: true
# Using case-insensitive name
sfcgalop -a "POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))" ISVALID
# Output: true
```

## Contributing

Contributions are welcome! Please follow these guidelines:

1. Use modern C++ features (C++17)
2. Follow existing code style (use clang-format)
3. Add tests for new operations
4. Update documentation

## License

This project is part of SFCGAL and follows the same license terms.

## Authors

SFCGAL Contributors

## See Also

- [SFCGAL Documentation](https://sfcgal.gitlab.io/SFCGAL/)
- [CGAL Documentation](https://www.cgal.org/)
- [OGC Simple Features](https://www.ogc.org/standards/sfa)
