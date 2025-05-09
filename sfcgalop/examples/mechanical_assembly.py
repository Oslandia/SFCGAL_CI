#!/usr/bin/env python3
# This file is meant to be executed via the sfcgal_dsl.py interpreter.
# Run with: python3 ../sfcgal_dsl.py <filename>
# The DSL automatically imports geometry constructors and operations.

# Mechanical example: assembly with revolution parts and transformations

# Main body - hollow cylinder (simulated by difference)
outer_cylinder = cylinder(height=20, radius=8)
inner_cylinder = cylinder(height=20, radius=6).translate([0, 0, 0])
main_body = outer_cylinder.difference(inner_cylinder)

# Fixing flanges at ends
flange_bottom = cylinder(height=2, radius=12).translate([0, 0, -1])
flange_top = cylinder(height=2, radius=12).translate([0, 0, 19])

# Fixing holes in flanges (simulated by small cylinders)
bolt_hole1 = cylinder(height=4, radius=1).translate([8, 0, -2])
bolt_hole2 = cylinder(height=4, radius=1).translate([-8, 0, -2])
bolt_hole3 = cylinder(height=4, radius=1).translate([0, 8, -2])
bolt_hole4 = cylinder(height=4, radius=1).translate([0, -8, -2])

# Same holes for upper flange
bolt_hole5 = cylinder(height=4, radius=1).translate([8, 0, 18])
bolt_hole6 = cylinder(height=4, radius=1).translate([-8, 0, 18])
bolt_hole7 = cylinder(height=4, radius=1).translate([0, 8, 18])
bolt_hole8 = cylinder(height=4, radius=1).translate([0, -8, 18])

# Decorative element - torus in the middle
decorative_ring = torus(major_radius=9, minor_radius=1).translate([0, 0, 10])

# Main assembly
assembly = main_body.union(flange_bottom).union(flange_top).union(decorative_ring)

# Subtraction of bolt holes
final_assembly = assembly.difference(bolt_hole1).difference(bolt_hole2).difference(bolt_hole3).difference(bolt_hole4).difference(bolt_hole5).difference(bolt_hole6).difference(bolt_hole7).difference(bolt_hole8)

# Rotation for presentation
final_assembly.rotate([15, 15, 0])