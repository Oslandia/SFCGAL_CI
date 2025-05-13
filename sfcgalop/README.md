# SFCGALOP - Command Line Utility for SFCGAL Operations

SFCGALOP is a command-line utility that provides access to the computational geometry operations in the SFCGAL library.

## Table of Contents

1. [Installation](#installation)
2. [Basic Usage](#basic-usage)
3. [Input Options](#input-options)
4. [Output Options](#output-options)
5. [Operations](#operations)
6. [Parameter Passing](#parameter-passing)
7. [Examples](#examples)
8. [Troubleshooting](#troubleshooting)

## Installation

SFCGALOP is built as part of the SFCGAL library. After building SFCGAL, the `sfcgalop` executable will be available in the build directory.

## Basic Usage

The basic syntax for running SFCGALOP is:

```
sfcgalop [OPTIONS] operation [operation_arg]
```

For example:

```bash
sfcgalop -a "POINT(0 0)" -b "POINT(1 1)" distance
```

This will calculate the distance between two points.

To see all available options:

```bash
sfcgalop --help
```

## Input Options

SFCGALOP supports various ways to specify input geometries:

### `-a, --geom-a=ARG`

Specifies the source for the first geometry (A). This is required for all operations. ARG can be:

- WKT string: `"POINT(0 0)"`
- WKB (hex): `"0101000000000000000000000000000000000000"`
- File path: `path/to/geometry.wkt` or `path/to/geometry.wkb`
- stdin: If set to `stdin`, SFCGALOP reads WKT from standard input
- stdin.wkb: If set to `stdin.wkb`, SFCGALOP reads WKB from standard input

If not specified, `stdin` is assumed.

### `-b, --geom-b=ARG`

Specifies the source for the second geometry (B). This is required for binary operations (operations that work on two geometries). The format options are the same as for `-a`.

## Output Options

### `-f, --format=ARG`

Specifies the output format. Available formats are:

- `wkt`: Well-Known Text (default)
- `wkb`: Well-Known Binary (hex encoded)
- `txt`: Plain text (currently same as WKT)
- `geojson`: GeoJSON format

Example:

```bash
sfcgalop -a "POINT(0 0)" -f geojson envelope
```

### `-p, --precision=N`

Sets the number of decimal places in output coordinates. Default is 6.

Example:

```bash
sfcgalop -a "POINT(0.123456789 0.987654321)" -p 3 centroid
```

This will output `POINT(0.123 0.988)`.

### `-q, --quiet`

Disables result output. Useful for benchmarking or when only interested in whether the operation succeeds.

### `-t, --time`

Print execution time in microseconds.

### `-v, --verbose`

Enable verbose output, showing more details about the operation being performed.

## Operations

To see all available operations:

```bash
sfcgalop --help
```

To get help on a specific operation:

```bash
sfcgalop --help-op=operation_name
```

For example:

```bash
sfcgalop --help-op=buffer3d
```

Operations are organized into categories:

- **Validity**: Operations to check geometry validity (`is_valid`, `is_validity_detail`, etc.)
- **Coordinate Handling**: Operations to manipulate coordinates (`drop_z`, `force_z`, etc.)
- **Spatial Analysis**: Spatial relationship operations (`intersects`, `covers`, etc.)
- **Construction**: Operations that create new geometries (`intersection`, `union`, etc.)
- **Metrics**: Operations that calculate measurements (`distance`, `area`, etc.)
- **Tesselation**: Operations for tesselating geometries (`tesselate`, `triangulate_2dz`, etc.)
- **Transformation**: Operations that transform geometries (`rotate`, `scale`, etc.)
- **Skeleton**: Operations related to skeletonization (`straight_skeleton`, etc.)
- **Alpha Shapes**: Operations related to alpha shapes
- **Visibility**: Operations related to visibility analysis

## Parameter Passing

Many operations accept additional parameters. There are two ways to specify parameters:

### 1. Positional Parameters

Parameters can be passed in a comma-separated list, in the order they are expected:

```bash
sfcgalop -a "LINESTRING(0 0, 1 1)" line_substring 0.25,0.75
```

This gets a substring of the linestring from 25% to 75% of its length.

### 2. Named Parameters

Parameters can also be passed by name:

```bash
sfcgalop -a "LINESTRING(0 0, 1 1)" line_substring start=0.25,end=0.75
```

This is equivalent to the previous example but uses named parameters.

Named parameters are especially useful when you want to specify only certain parameters or specify them in a different order:

```bash
sfcgalop -a "POLYGON((0 0, 0 10, 10 10, 10 0, 0 0))" offset_polygon radius=5.0
```

## Examples

Here are some common usage examples:

### Calculate the distance between two points:

```bash
sfcgalop -a "POINT(0 0)" -b "POINT(3 4)" distance
```

Output: `5`

### Create a buffer around a point:

```bash
sfcgalop -a "POINT(0 0)" -f wkt buffer3d radius=2.0,segments=32,buffer_type=0
```

### Check if a polygon is valid:

```bash
sfcgalop -a "POLYGON((0 0, 0 10, 10 10, 10 0, 0 0))" is_valid
```

Output: `true`

### Get detailed validity information:

```bash
sfcgalop -a "POLYGON((0 0, 10 0, 10 10, 0 10, 5 5, 0 0))" is_validity_detail
```

### Process geometries from files:

```bash
sfcgalop -a input.wkt -f geojson convexhull
```

### Pipe WKT data:

```bash
echo "POINT(0 0)" | sfcgalop -f wkt buffer3d radius=1.0,segments=16,buffer_type=0
```

### Rotate a geometry:

```bash
sfcgalop -a "POLYGON((0 0, 0 1, 1 1, 1 0, 0 0))" rotate angle=1.5708
```

This rotates the polygon 90 degrees (π/2 radians).

### Scale a geometry:

```bash
sfcgalop -a "POLYGON((0 0, 0 1, 1 1, 1 0, 0 0))" scale_3d sx=2.0,sy=3.0,sz=1.0
```

This scales the polygon by 2x in X direction and 3x in Y direction.

## Troubleshooting

### Common Errors

1. **"Operation requires a second geometry"**: This error occurs when you try to run a binary operation without specifying the second geometry (`-b`). Make sure to provide both geometries for operations that require two inputs.

2. **"Invalid parameter value(s)"**: This occurs when the parameters passed to an operation are not valid. Check the parameter names and values.

3. **"Operation failed"**: This generic error might be followed by a more specific message. It indicates that the operation could not be completed successfully.

4. **Parsing errors**: If your WKT or WKB input cannot be parsed, you'll get an error. Make sure your input format is correct.

### Tips

- Use the `-v` (verbose) flag to get more information about what's happening
- Check operation help using `--help-op=operation_name` to understand required parameters
- For binary operations, make sure the order of geometries is correct, as some operations are not commutative (e.g., `difference`)
- If working with large geometries, consider storing them in files rather than passing them directly on the command line
