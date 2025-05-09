#!/usr/bin/env python3
# This file is meant to be executed via the sfcgal_dsl.py interpreter.
# Run with: python3 ../sfcgal_dsl.py <filename>
# The DSL automatically imports geometry constructors and operations.

# For standalone execution, uncomment this:
# from sfcgal_dsl import point, linestring, polygon

# Example: basic 2D geometries
# Simple point
p = point([10, 20])

# 3D LineString
ls = linestring([[0, 0, 0], [10, 0, 5], [20, 10, 8]])

# Polygon with hole
poly = polygon([
    [[0, 0], [20, 0], [20, 20], [0, 20], [0, 0]],  # Exterior
    [[5, 5], [15, 5], [15, 15], [5, 15], [5, 15], [5, 5]]   # Hole
])

# Return the polygon
poly