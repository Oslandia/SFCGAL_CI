#!/bin/bash

# SFCGALOP Usage Examples
# This file contains various examples showing how to use sfcgalop

echo "=== SFCGALOP Examples ==="
echo

# Basic area calculation
echo "1. Calculate polygon area:"
echo "   sfcgalop -a \"POLYGON((0 0, 0 10, 10 10, 10 0, 0 0))\" area"
./build/sfcgalop/sfcgalop -a "POLYGON((0 0, 0 10, 10 10, 10 0, 0 0))" area
echo

# Distance between points
echo "2. Distance between two points:"
echo "   sfcgalop -a \"POINT(0 0)\" -b \"POINT(3 4)\" distance"
./build/sfcgalop/sfcgalop -a "POINT(0 0)" -b "POINT(3 4)" distance
echo

# Intersection of polygons
echo "3. Intersection of two polygons:"
echo "   sfcgalop -a \"POLYGON((0 0, 4 0, 4 4, 0 4, 0 0))\" -b \"POLYGON((2 2, 6 2, 6 6, 2 6, 2 2))\" intersection"
./build/sfcgalop/sfcgalop -a "POLYGON((0 0, 4 0, 4 4, 0 4, 0 0))" -b "POLYGON((2 2, 6 2, 6 6, 2 6, 2 2))" intersection
echo

# Convex hull
echo "4. Convex hull of multipoint:"
echo "   sfcgalop -a \"MULTIPOINT((0 0),(1 1),(1 0),(0 1),(0.5 0.5))\" convexhull"
./build/sfcgalop/sfcgalop -a "MULTIPOINT((0 0),(1 1),(1 0),(0 1),(0.5 0.5))" convexhull
echo

# Validation example
echo "5. Validate geometry (invalid polygon):"
echo "   sfcgalop --validate -a \"POLYGON((0 0, 0 10, 10 10, 10 0, 1 0))\" area"
./build/sfcgalop/sfcgalop --validate -a "POLYGON((0 0, 0 10, 10 10, 10 0, 1 0))" area 2>&1 | head -10
echo

# 3D operations
echo "6. Volume of a cube:"
echo "   sfcgalop -a \"SOLID((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),((0 0 0,0 0 1,0 1 1,0 1 0,0 0 0)),((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0)),((1 1 1,0 1 1,0 0 1,1 0 1,1 1 1)),((1 1 1,1 0 1,1 0 0,1 1 0,1 1 1)),((1 1 1,1 1 0,0 1 0,0 1 1,1 1 1))))\" volume"
./build/sfcgalop/sfcgalop -a "SOLID((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),((0 0 0,0 0 1,0 1 1,0 1 0,0 0 0)),((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0)),((1 1 1,0 1 1,0 0 1,1 0 1,1 1 1)),((1 1 1,1 0 1,1 0 0,1 1 0,1 1 1)),((1 1 1,1 1 0,0 1 0,0 1 1,1 1 1))))" volume
echo

# Transformation - translate
echo "7. Translate a point:"
echo "   sfcgalop -a \"POINT(0 0)\" translate 10,20"
./build/sfcgalop/sfcgalop -a "POINT(0 0)" translate 10,20
echo

# Simplify geometry
echo "8. Simplify linestring:"
echo "   sfcgalop -a \"LINESTRING(0 0,0.1 0.05,0.2 0.1,0.3 0.15,1 1)\" simplify 0.2"
./build/sfcgalop/sfcgalop -a "LINESTRING(0 0,0.1 0.05,0.2 0.1,0.3 0.15,1 1)" simplify 0.2
echo

# Buffer operation
echo "9. 3D buffer around a point:"
echo "   sfcgalop -a \"POINT(0 0 0)\" buffer3d 2.0"
./build/sfcgalop/sfcgalop -a "POINT(0 0 0)" buffer3d 2.0 | head -1
echo "[Truncated for brevity - generates a sphere]"
echo

# Using stdin
echo "10. Reading from stdin:"
echo "    echo \"POLYGON((0 0, 0 1, 1 1, 1 0, 0 0))\" | sfcgalop -a stdin area"
echo "POLYGON((0 0, 0 1, 1 1, 1 0, 0 0))" | ./build/sfcgalop/sfcgalop -a stdin area
echo

# GeoJSON output
echo "11. Output as GeoJSON:"
echo "    sfcgalop -a \"POINT(1 2)\" --format geojson centroid"
./build/sfcgalop/sfcgalop -a "POINT(1 2)" --format geojson centroid
echo

# Timing operations
echo "12. Show execution time:"
echo "    sfcgalop -t -a \"POLYGON((0 0, 0 1000, 1000 1000, 1000 0, 0 0))\" area"
./build/sfcgalop/sfcgalop -t -a "POLYGON((0 0, 0 1000, 1000 1000, 1000 0, 0 0))" area
echo

echo "=== Advanced Examples ==="
echo

# Straight skeleton
echo "13. Straight skeleton of a polygon:"
echo "    sfcgalop -a \"POLYGON((0 0, 4 0, 4 4, 2 4, 2 2, 0 2, 0 0))\" straightskeleton"
./build/sfcgalop/sfcgalop -a "POLYGON((0 0, 4 0, 4 4, 2 4, 2 2, 0 2, 0 0))" straightskeleton | head -1
echo "[Complex output truncated]"
echo

# Minkowski sum
echo "14. Minkowski sum:"
echo "    sfcgalop -a \"POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))\" -b \"POLYGON((0 0, 0.5 0, 0.5 0.5, 0 0.5, 0 0))\" minkowskisum"
./build/sfcgalop/sfcgalop -a "POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))" -b "POLYGON((0 0, 0.5 0, 0.5 0.5, 0 0.5, 0 0))" minkowskisum
echo

# Collection operations
echo "15. Extract polygons from collection:"
echo "    sfcgalop -a \"GEOMETRYCOLLECTION(POINT(0 0),POLYGON((0 0,1 0,1 1,0 1,0 0)),LINESTRING(2 2,3 3))\" collection_extract 3"
./build/sfcgalop/sfcgalop -a "GEOMETRYCOLLECTION(POINT(0 0),POLYGON((0 0,1 0,1 1,0 1,0 0)),LINESTRING(2 2,3 3))" collection_extract 3
echo

echo "=== List all operations ==="
echo "To see all available operations, run:"
echo "   sfcgalop --list"
echo
echo "For help on a specific operation, run:"
echo "   sfcgalop --help-op=<operation>"