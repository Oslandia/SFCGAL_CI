#!/usr/bin/env python3
"""
SFCGAL DSL - OpenSCAD-style Domain Specific Language for SFCGAL
Enables writing geometric scripts in a fluent style similar to OpenSCAD
"""

import re
import subprocess
import sys
import os
import tempfile
from typing import List, Dict, Any, Optional, Union
import json

class GeometryObject:
    """Represents a geometric object with fluent methods"""

    def __init__(self, geometry_data: str = "", operation_type: str = "geometry"):
        self.geometry_data = geometry_data
        self.operation_type = operation_type
        self.sfcgalop_path = "./build/sfcgalop/sfcgalop"

    def translate(self, coords: Union[List[float], Dict[str, float]]) -> 'GeometryObject':
        """Apply a translation to the geometry"""
        if isinstance(coords, list):
            if len(coords) >= 3:
                x, y, z = coords[0], coords[1], coords[2]
            elif len(coords) == 2:
                x, y, z = coords[0], coords[1], 0.0
            else:
                raise ValueError("translate() requires at least 2 coordinates")
        else:
            x = coords.get('x', 0.0)
            y = coords.get('y', 0.0)
            z = coords.get('z', 0.0)

        # Execute the transformation
        result_data = self._execute_transformation("translate", f"x={x},y={y},z={z}")
        return GeometryObject(result_data, "transformed")

    def rotate(self, coords: Union[List[float], Dict[str, float]]) -> 'GeometryObject':
        """Apply a rotation to the geometry (in degrees)"""
        if isinstance(coords, list):
            if len(coords) >= 3:
                x, y, z = coords[0], coords[1], coords[2]
            elif len(coords) == 2:
                x, y, z = coords[0], coords[1], 0.0
            else:
                x, y, z = coords[0], 0.0, 0.0
        else:
            x = coords.get('x', 0.0)
            y = coords.get('y', 0.0)
            z = coords.get('z', 0.0)

        # Currently only support Z-axis rotation
        if x != 0.0 or y != 0.0:
            print("Warning: Only Z-axis rotation currently supported", file=sys.stderr)

        result_data = self._execute_transformation("rotate", f"angle={z},axis=z")
        return GeometryObject(result_data, "transformed")

    def scale(self, factor: Union[float, List[float], Dict[str, float]]) -> 'GeometryObject':
        """Apply scaling to the geometry"""
        if isinstance(factor, (int, float)):
            params = f"s={factor}"
        elif isinstance(factor, list):
            if len(factor) >= 3:
                params = f"sx={factor[0]},sy={factor[1]},sz={factor[2]}"
            elif len(factor) == 2:
                params = f"sx={factor[0]},sy={factor[1]},sz=1.0"
            else:
                params = f"s={factor[0]}"
        else:
            sx = factor.get('sx', factor.get('factor', 1.0))
            sy = factor.get('sy', factor.get('factor', 1.0))
            sz = factor.get('sz', factor.get('factor', 1.0))
            params = f"sx={sx},sy={sy},sz={sz}"

        result_data = self._execute_transformation("scale", params)
        return GeometryObject(result_data, "transformed")

    def union(self, other: 'GeometryObject') -> 'GeometryObject':
        """Boolean union with another geometry"""
        result_data = self._execute_boolean_operation("union3d", other)
        return GeometryObject(result_data, "boolean")

    def intersection(self, other: 'GeometryObject') -> 'GeometryObject':
        """Boolean intersection with another geometry"""
        result_data = self._execute_boolean_operation("intersection3d", other)
        return GeometryObject(result_data, "boolean")

    def difference(self, other: 'GeometryObject') -> 'GeometryObject':
        """Boolean difference with another geometry"""
        result_data = self._execute_boolean_operation("difference3d", other)
        return GeometryObject(result_data, "boolean")

    def to_solid(self) -> 'GeometryObject':
        """Convert a PolyhedralSurface to a Solid"""
        result_data = self._execute_transformation("to_solid", "")
        return GeometryObject(result_data, "solid")

    def _execute_transformation(self, operation: str, params: str) -> str:
        """Execute a transformation via sfcgalop"""
        try:
            cmd = [self.sfcgalop_path, "-a", "stdin", operation, params]
            process = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                     text=True)
            result_stdout, result_stderr = process.communicate(input=self.geometry_data, timeout=30)

            if process.returncode == 0:
                return result_stdout.strip()
            else:
                raise RuntimeError(f"SFCGAL {operation} error: {result_stderr}")
        except subprocess.TimeoutExpired:
            raise RuntimeError(f"Timeout during {operation} operation")
        except Exception as e:
            raise RuntimeError(f"Error during {operation}: {str(e)}")

    def _execute_boolean_operation(self, operation: str, other: 'GeometryObject') -> str:
        """Execute a boolean operation via sfcgalop"""
        f1_name = None
        f2_name = None
        try:
            # Write geometries to temporary files
            # Close handles before spawning subprocess to avoid Windows sharing violations
            with tempfile.NamedTemporaryFile(mode='w', suffix='.wkt', delete=False) as f1:
                f1.write(self.geometry_data)
                f1.flush()
                f1_name = f1.name

            with tempfile.NamedTemporaryFile(mode='w', suffix='.wkt', delete=False) as f2:
                f2.write(other.geometry_data)
                f2.flush()
                f2_name = f2.name

            # Execute subprocess with files closed (no sharing violations on Windows)
            cmd = [self.sfcgalop_path, "-a", f1_name, "-b", f2_name, operation]
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            result_stdout, result_stderr = process.communicate(timeout=30)

            if process.returncode == 0:
                return result_stdout.strip()
            else:
                # Boolean operations may not be supported
                print(f"Warning: {operation} operation may not be supported yet: {result_stderr}", file=sys.stderr)
                return self.geometry_data  # Return original geometry

        except Exception as e:
            print(f"Warning: Boolean operation {operation} failed: {str(e)}", file=sys.stderr)
            return self.geometry_data  # Return original geometry
        finally:
            # Cleanup temporary files in finally block to ensure removal
            if f1_name:
                try:
                    os.unlink(f1_name)
                except PermissionError:
                    print(f"Warning: Could not remove temporary file {f1_name}", file=sys.stderr)
                except OSError:
                    pass  # File may already be deleted
            if f2_name:
                try:
                    os.unlink(f2_name)
                except PermissionError:
                    print(f"Warning: Could not remove temporary file {f2_name}", file=sys.stderr)
                except OSError:
                    pass  # File may already be deleted

    def to_wkt(self) -> str:
        """Return the WKT representation of the geometry"""
        return self.geometry_data

    def __str__(self) -> str:
        return self.to_wkt()

# Primitive geometry creation functions
def sphere(radius: float = 1.0, **kwargs) -> GeometryObject:
    """Create a sphere"""
    # Parameter validation
    if radius <= 0:
        raise ValueError("sphere() radius must be positive")

    params = [f"radius={radius}"]

    # Optional parameters with validation
    if 'x' in kwargs:
        params.append(f"x={float(kwargs['x'])}")
    if 'y' in kwargs:
        params.append(f"y={float(kwargs['y'])}")
    if 'z' in kwargs:
        params.append(f"z={float(kwargs['z'])}")
    if 'num_vertical' in kwargs:
        num_vertical = int(kwargs['num_vertical'])
        if num_vertical < 3:
            raise ValueError("sphere() num_vertical must be at least 3")
        params.append(f"num_vertical={num_vertical}")
    if 'num_horizontal' in kwargs:
        num_horizontal = int(kwargs['num_horizontal'])
        if num_horizontal < 3:
            raise ValueError("sphere() num_horizontal must be at least 3")
        params.append(f"num_horizontal={num_horizontal}")

    return _create_primitive("make_sphere", ",".join(params))

def cube(size: Union[float, List[float]] = 1.0) -> GeometryObject:
    """Create a cube or box"""
    if isinstance(size, (int, float)):
        if size <= 0:
            raise ValueError("cube() size must be positive")
        return _create_primitive("make_cube", f"size={size}")
    elif isinstance(size, list) and len(size) >= 3:
        if any(s <= 0 for s in size[:3]):
            raise ValueError("cube() dimensions must be positive")
        return _create_primitive("make_box", f"x_extent={size[0]},y_extent={size[1]},z_extent={size[2]}")
    else:
        raise ValueError("cube() requires a single size value or a list of [x, y, z] dimensions")

def box(dimensions: List[float]) -> GeometryObject:
    """Create a box with specific dimensions"""
    if len(dimensions) < 3:
        raise ValueError("box() requires [x, y, z] dimensions")
    if any(d <= 0 for d in dimensions[:3]):
        raise ValueError("box() dimensions must be positive")
    return _create_primitive("make_box", f"x_extent={dimensions[0]},y_extent={dimensions[1]},z_extent={dimensions[2]}")

def cylinder(height: float = 1.0, radius: float = 1.0, **kwargs) -> GeometryObject:
    """Create a cylinder"""
    # Parameter validation
    if height <= 0:
        raise ValueError("cylinder() height must be positive")
    if radius <= 0:
        raise ValueError("cylinder() radius must be positive")

    params = [f"height={height}", f"radius={radius}"]

    if 'num_radial' in kwargs:
        num_radial = int(kwargs['num_radial'])
        if num_radial < 3:
            raise ValueError("cylinder() num_radial must be at least 3")
        params.append(f"num_radial={num_radial}")

    return _create_primitive("make_cylinder", ",".join(params))

def cone(height: float = 1.0, bottom_radius: float = None, top_radius: float = 0.0,
         radius: float = None, **kwargs) -> GeometryObject:
    """Create a cone (simple or truncated)"""
    # Support both syntaxes: cone(radius=...) and cone(bottom_radius=...)
    if bottom_radius is None:
        bottom_radius = radius if radius is not None else 1.0

    # Parameter validation
    if height <= 0:
        raise ValueError("cone() height must be positive")
    if bottom_radius < 0 or top_radius < 0:
        raise ValueError("cone() radii must be non-negative")
    if bottom_radius == 0 and top_radius == 0:
        raise ValueError("cone() requires at least one non-zero radius")

    params = [f"height={height}", f"bottom_radius={bottom_radius}", f"top_radius={top_radius}"]

    if 'num_radial' in kwargs:
        num_radial = int(kwargs['num_radial'])
        if num_radial < 3:
            raise ValueError("cone() num_radial must be at least 3")
        params.append(f"num_radial={num_radial}")

    return _create_primitive("make_cone", ",".join(params))

def torus(major_radius: float = 2.0, minor_radius: float = 0.5, **kwargs) -> GeometryObject:
    """Create a torus"""
    # Parameter validation
    if major_radius <= 0:
        raise ValueError("torus() major_radius must be positive")
    if minor_radius <= 0:
        raise ValueError("torus() minor_radius must be positive")
    if minor_radius >= major_radius:
        raise ValueError("torus() minor_radius must be less than major_radius")

    params = [f"major_radius={major_radius}", f"minor_radius={minor_radius}"]

    if 'num_major' in kwargs:
        num_major = int(kwargs['num_major'])
        if num_major < 3:
            raise ValueError("torus() num_major must be at least 3")
        params.append(f"num_major={num_major}")
    if 'num_minor' in kwargs:
        num_minor = int(kwargs['num_minor'])
        if num_minor < 3:
            raise ValueError("torus() num_minor must be at least 3")
        params.append(f"num_minor={num_minor}")

    return _create_primitive("make_torus", ",".join(params))

# Basic geometry functions
def point(coords: List[float]) -> GeometryObject:
    """Create a point"""
    if not isinstance(coords, list):
        raise ValueError("point() requires a list of coordinates")
    if len(coords) == 2:
        wkt = f"POINT({coords[0]} {coords[1]})"
    elif len(coords) == 3:
        wkt = f"POINT Z({coords[0]} {coords[1]} {coords[2]})"
    else:
        raise ValueError("point() requires 2 or 3 coordinates")
    return GeometryObject(wkt, "basic")

def linestring(coords: List[List[float]]) -> GeometryObject:
    """Create a LineString"""
    if not isinstance(coords, list) or len(coords) < 2:
        raise ValueError("linestring() requires at least 2 points")

    if not all(isinstance(pt, list) for pt in coords):
        raise ValueError("linestring() coordinates must be lists of numbers")

    # Check dimension consistency
    first_dim = len(coords[0]) if coords else 0
    if not all(len(pt) == first_dim for pt in coords):
        raise ValueError("linestring() all points must have the same dimension")

    if first_dim == 2:
        coord_strs = [f"{pt[0]} {pt[1]}" for pt in coords]
        wkt = f"LINESTRING({','.join(coord_strs)})"
    elif first_dim == 3:
        coord_strs = [f"{pt[0]} {pt[1]} {pt[2]}" for pt in coords]
        wkt = f"LINESTRING Z({','.join(coord_strs)})"
    else:
        raise ValueError("linestring() coordinates must be 2D or 3D")

    return GeometryObject(wkt, "basic")

def polygon(coords: Union[List[List[float]], List[List[List[float]]]]) -> GeometryObject:
    """Create a polygon (with optional holes)"""
    if not coords:
        raise ValueError("polygon() requires coordinates")

    # Determine if it's a polygon with holes
    if isinstance(coords[0][0], list):
        # Format: [exterior_ring, hole1, hole2, ...]
        rings = coords
    else:
        # Simple format: [point1, point2, ...]
        rings = [coords]

    # Ring validation
    for i, ring in enumerate(rings):
        if len(ring) < 3:
            ring_type = "exterior ring" if i == 0 else f"hole {i}"
            raise ValueError(f"polygon() {ring_type} requires at least 3 points")

        # Check that ring is closed
        if ring[0] != ring[-1]:
            ring_type = "exterior ring" if i == 0 else f"hole {i}"
            raise ValueError(f"polygon() {ring_type} must be closed (first and last points must be equal)")

    # Build WKT
    if len(rings[0][0]) == 2:
        ring_strs = []
        for ring in rings:
            coord_strs = [f"{pt[0]} {pt[1]}" for pt in ring]
            ring_strs.append(f"({','.join(coord_strs)})")
        wkt = f"POLYGON({','.join(ring_strs)})"
    else:
        ring_strs = []
        for ring in rings:
            coord_strs = [f"{pt[0]} {pt[1]} {pt[2]}" for pt in ring]
            ring_strs.append(f"({','.join(coord_strs)})")
        wkt = f"POLYGON Z({','.join(ring_strs)})"

    return GeometryObject(wkt, "basic")

def _create_primitive(command: str, params: str) -> GeometryObject:
    """Helper function to create primitives"""
    try:
        sfcgalop_path = "./build/sfcgalop/sfcgalop"
        cmd = [sfcgalop_path, command, params]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)

        if result.returncode == 0:
            return GeometryObject(result.stdout.strip(), "primitive")
        else:
            raise RuntimeError(f"SFCGAL {command} error: {result.stderr}")
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"Timeout during {command} creation")
    except Exception as e:
        raise RuntimeError(f"Error creating {command}: {str(e)}")

class SFCGALDSLParser:
    """Simple parser for Python-style DSL expressions"""

    def __init__(self):
        # Expose geometry functions in global namespace
        self.globals = {
            'sphere': sphere,
            'cube': cube,
            'box': box,
            'cylinder': cylinder,
            'cone': cone,
            'torus': torus,
            'point': point,
            'linestring': linestring,
            'polygon': polygon,
        }

    def parse_and_execute(self, code: str) -> GeometryObject:
        """Parse and execute DSL code"""
        try:
            # Remove comments
            code = re.sub(r'//.*', '', code)
            code = re.sub(r'/\*.*?\*/', '', code, flags=re.DOTALL)

            # Execute Python code
            local_vars = {}
            exec(code, self.globals, local_vars)

            # Find the last expression (result)
            lines = [line.strip() for line in code.strip().split('\n') if line.strip()]
            if lines:
                last_line = lines[-1]
                if not last_line.endswith(';'):
                    # Evaluate the last expression
                    result = eval(last_line, self.globals, local_vars)
                    if isinstance(result, GeometryObject):
                        return result

            # If no return found, search in local variables
            for var_name, var_value in local_vars.items():
                if isinstance(var_value, GeometryObject):
                    return var_value

            raise ValueError("No geometry result found in DSL code")

        except Exception as e:
            raise RuntimeError(f"DSL parsing/execution error: {str(e)}")

def main():
    """Main entry point"""
    if len(sys.argv) != 2:
        print("Usage: python3 sfcgal_dsl.py <script.py>")
        print("\nExamples:")
        print("  # Simple sphere")
        print("  sphere(radius=10)")
        print()
        print("  # Transformed geometry")
        print("  sphere(radius=5).translate([10, 0, 0]).scale(2)")
        print()
        print("  # Boolean operations")
        print("  sphere(radius=10).union(cube(size=8).translate([15, 0, 0]))")
        sys.exit(1)

    script_file = sys.argv[1]

    if not os.path.exists(script_file):
        print(f"Error: File {script_file} not found")
        sys.exit(1)

    try:
        # Read the script
        with open(script_file, 'r') as f:
            code = f.read()

        # Parse and execute the code
        parser = SFCGALDSLParser()
        result = parser.parse_and_execute(code)

        # Display WKT result
        print(result.to_wkt())

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
