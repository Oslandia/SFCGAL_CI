#!/usr/bin/env python3
# Simple example: sphere with transformations
# This file is meant to be executed via the sfcgal_dsl.py interpreter.
# Run with: python3 sfcgal_dsl.py examples/simple_sphere.py
# The DSL automatically imports and provides: sphere, translate, scale, etc.

sphere(radius=10).translate([5, 0, 2]).scale(1.5)