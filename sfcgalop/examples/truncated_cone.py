#!/usr/bin/env python3
# This file is meant to be executed via the sfcgal_dsl.py interpreter.
# Run with: python3 ../sfcgal_dsl.py <filename>
# The DSL automatically imports geometry constructors and operations.

# Example: truncated cone (frustum)
cone(height=8, bottom_radius=5, top_radius=2, num_radial=12)