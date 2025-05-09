#!/usr/bin/env python3
# This file is meant to be executed via the sfcgal_dsl.py interpreter.
# Run with: python3 ../sfcgal_dsl.py <filename>
# The DSL automatically imports geometry constructors and operations.

# Complex example: architectural structure with boolean operations and transformations

# Main base - slab
base = cube([40, 30, 3])

# Main tower at center
main_tower = cylinder(height=25, radius=4).translate([20, 15, 3])

# Smaller corner towers
corner_tower1 = cylinder(height=15, radius=2).translate([5, 5, 3])
corner_tower2 = cylinder(height=15, radius=2).translate([35, 5, 3])
corner_tower3 = cylinder(height=15, radius=2).translate([5, 25, 3])
corner_tower4 = cylinder(height=15, radius=2).translate([35, 25, 3])

# Connections between towers (walls)
wall1 = cube([30, 2, 12]).translate([5, 14, 3])
wall2 = cube([2, 20, 12]).translate([19, 5, 3])

# Truncated cone roof on main tower
roof = cone(height=8, bottom_radius=5, top_radius=1).translate([20, 15, 28])

# Cutouts for doors and windows in walls
door_opening = cube([3, 1, 8]).translate([18, 14, 3])
window1 = cube([2, 1, 3]).translate([10, 14, 8])
window2 = cube([2, 1, 3]).translate([25, 14, 8])

# Assembly with boolean operations
# Union of all solid elements
structure = base.union(main_tower).union(corner_tower1).union(corner_tower2).union(corner_tower3).union(corner_tower4).union(wall1).union(wall2).union(roof)

# Subtraction of openings
final_structure = structure.difference(door_opening).difference(window1).difference(window2)

# Return final structure
final_structure