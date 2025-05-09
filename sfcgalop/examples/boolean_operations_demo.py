#!/usr/bin/env python3
# This file is meant to be executed via the sfcgal_dsl.py interpreter.
# Run with: python3 ../sfcgal_dsl.py <filename>
# The DSL automatically imports geometry constructors and operations.

# For standalone execution, uncomment this:
# from sfcgal_dsl import sphere, cube, cylinder, torus

# Demonstration of boolean operations with primitive geometries

# Create two overlapping spheres
sphere1 = sphere(radius=8, x=0, y=0, z=0)
sphere2 = sphere(radius=8, x=10, y=0, z=0)

# Union: merge the two spheres
union_result = sphere1.union(sphere2)

# Intersection: common volume
intersection_result = sphere1.intersection(sphere2)

# Difference: first sphere minus the second
difference_result = sphere1.difference(sphere2)

# Example with cube and cylinder
base_cube = cube([20, 20, 10])
cutting_cylinder = cylinder(height=12, radius=6).translate([10, 10, -1])

# Drilling the cube with the cylinder
perforated_cube = base_cube.difference(cutting_cylinder)

# Complex combination: torus intersecting with a box
torus_shape = torus(major_radius=12, minor_radius=3)
bounding_box = cube([20, 20, 8]).translate([-10, -10, -4])
constrained_torus = torus_shape.intersection(bounding_box)

# Final transformation with scaling and rotation
final_shape = constrained_torus.scale([1.2, 1.2, 0.8]).rotate([0, 0, 45])

# Return final shape
final_shape