#!/usr/bin/env python3
# This file is meant to be executed via the sfcgal_dsl.py interpreter.
# Run with: python3 ../sfcgal_dsl.py <filename>
# The DSL automatically imports geometry constructors and operations.

# For standalone execution, uncomment this:
# from sfcgal_dsl import point, linestring, polygon, cube, cylinder, cone

# Mixed example: basic geometries and volumes with transformations

# Basic 2D geometries
central_point = point([0, 0])
reference_line = linestring([[-10, -5], [0, 0], [10, 5], [15, 10]])

# Complex polygon with hole
outer_boundary = [
    [-15, -10], [15, -10], [15, 10], [-15, 10], [-15, -10]
]
inner_hole = [
    [-5, -3], [5, -3], [5, 3], [-5, 3], [-5, -3]
]
complex_polygon = polygon([outer_boundary, inner_hole])

# 3D volumes based conceptually on 2D shapes
# Simulated extrusion of polygon via box with cutout
base_volume = cube([30, 20, 5]).translate([-15, -10, 0])
hole_volume = cube([10, 6, 6]).translate([-5, -3, -1])
extruded_shape = base_volume.difference(hole_volume)

# Add cylindrical elements at corners
corner_posts = [
    cylinder(height=8, radius=1).translate([-15, -10, 5]),
    cylinder(height=8, radius=1).translate([15, -10, 5]),
    cylinder(height=8, radius=1).translate([15, 10, 5]),
    cylinder(height=8, radius=1).translate([-15, 10, 5])
]

# Assembly with posts
structure_with_posts = extruded_shape
for post in corner_posts:
    structure_with_posts = structure_with_posts.union(post)

# Add central decorative element
central_feature = cone(
    height=6,
    bottom_radius=3,
    top_radius=1,
    num_radial=8
).translate([0, 0, 5])

# Final result with transformation
final_result = structure_with_posts.union(central_feature).rotate([0, 0, 30]).scale([1.1, 1.1, 1.0])

# Return final result
final_result