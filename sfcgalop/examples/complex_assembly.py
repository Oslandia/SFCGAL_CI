#!/usr/bin/env python3
# This file is meant to be executed via the sfcgal_dsl.py interpreter.
# Run with: python3 ../sfcgal_dsl.py <filename>
# The DSL automatically imports geometry constructors and operations.

# For standalone execution, uncomment this:
# from sfcgal_dsl import cube, cylinder

# Complex example: assembly with multiple transformations

# Large base
base = cube([30, 30, 5])

# Corner pillars
pillar1 = cylinder(height=15, radius=3).translate([5, 5, 5])
pillar2 = cylinder(height=15, radius=3).translate([25, 5, 5])
pillar3 = cylinder(height=15, radius=3).translate([5, 25, 5])
pillar4 = cylinder(height=15, radius=3).translate([25, 25, 5])

# Roof
roof = cube([30, 30, 3]).translate([0, 0, 20])

# Final assembly (for now just return the base)
base