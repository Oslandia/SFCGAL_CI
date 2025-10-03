# SFCGAL WebAssembly Binding

A complete WebAssembly binding for SFCGAL (Simple Features for CGAL), enabling 3D geometric operations in web browsers using the actual SFCGAL C++ library compiled to WebAssembly.

## Features

- **Real SFCGAL Library**: Uses the actual SFCGAL C++ implementation, not a simplified JavaScript port
- **Full 3D Support**: All SFCGAL geometric operations including 3D primitives (Sphere, Cylinder, Cone, Torus)
- **Comprehensive Operations**: Distance, intersection, union, difference, buffer, extrusion, and more
- **Modern Web Integration**: ES6 modules, TypeScript-friendly, runs in browsers and Node.js

## Quick Start

### Prerequisites

- Linux or macOS
- CMake >= 3.10
- C++17 compatible compiler
- Git, curl, tar, xz

### Build Instructions

1. **Install Emscripten** (if not already installed):
```bash
./scripts/install-emscripten.sh
```

2. **Build dependencies** (GMP, MPFR, Boost, CGAL):

   ```bash
   ./scripts/build-deps.sh
   ```
   Builds dependencies with a cmake files to avoid modifying CGAL files.

3. **Build SFCGAL WebAssembly library**:

   ```bash
   ./scripts/build-sfcgal-wasm.sh
   ```


4. **Build the binding**:
```bash
./scripts/build.sh
```
This script automatically detects which SFCGAL library is available and uses it to build the final WebAssembly module (takes ~20 seconds).

> **Note:** See `cmake/README.md` for detailed comparison of both build methods.

### Quick Build (All-in-One)

```bash
yes | ./scripts/build-complete.sh
```

### Output

The build produces two files in `examples/`:
- `sfcgal.js` - JavaScript loader and wrapper (~107KB)
- `sfcgal.wasm` - WebAssembly binary module (~3.0MB)

### Run Examples

Start a local web server:
```bash
./scripts/serve.sh
```

Then open http://localhost:8000 in your browser to see the examples.

## Project Structure

```
sfcgal-wasm/
├── scripts/
│   ├── install-emscripten.sh    # Install Emscripten SDK
│   ├── build-deps.sh             # Build WebAssembly dependencies
│   ├── build-sfcgal-wasm.sh      # Build SFCGAL static library
│   ├── build.sh                  # Build the WebAssembly binding
│   └── serve.sh                  # Start development server
├── src/
│   └── sfcgal-binding.cpp        # C++ binding using Embind
├── examples/
│   ├── index.html                # Examples index with organized categories
│   ├── sfcgal-module.js          # Shared SFCGAL module loader
│   ├── canvas2d-renderer.js      # 2D Canvas renderer utility
│   ├── three3d-renderer.js       # 3D Three.js renderer utility
│   │
│   ├── 2D Examples (Canvas 2D):
│   │   ├── area-2d.html          # Area calculation
│   │   ├── distance-2d.html      # Distance between geometries
│   │   ├── centroid-2d.html      # Centroid calculation
│   │   ├── intersection-2d.html  # Intersection operations
│   │   ├── union-2d.html         # Union operations
│   │   ├── difference-2d.html    # Difference operations
│   │   ├── convexhull-2d.html    # Convex hull
│   │   ├── buffer-2d.html        # Buffer operations
│   │   └── transformations-2d.html # Translate, rotate, scale
│   │
│   ├── 3D Examples (Three.js):
│   │   ├── extrude-3d.html       # 3D extrusion with wireframe mode
│   │   ├── area-3d.html          # 3D area calculation
│   │   ├── distance-3d.html      # 3D distance
│   │   ├── boolean-3d.html       # 3D boolean ops (intersection3D, union3D, difference3D)
│   │   └── transformations-3d.html # 3D transformations
│   │
└── deps/                         # Dependencies (built by scripts)
```

## API Usage

### Basic Example

```javascript
// Using the shared module loader (recommended)
import { loadSFCGAL } from './sfcgal-module.js';

const sfcgal = await loadSFCGAL();

// Calculate distance between two geometries
const wkt1 = "POINT(0 0)";
const wkt2 = "POINT(3 4)";
const distance = sfcgal.distance(wkt1, wkt2);
console.log(`Distance: ${distance}`); // 5.0

// Union of polygons
const poly1 = "POLYGON((0 0,10 0,10 10,0 10,0 0))";
const poly2 = "POLYGON((5 5,15 5,15 15,5 15,5 5))";
const result = sfcgal.union(poly1, poly2);

// Extrude to 3D
const polygon = "POLYGON((0 0,1 0,1 1,0 1,0 0))";
const solid = sfcgal.extrude(polygon, 0, 0, 2);
```

### Direct Module Usage

```javascript
import SFCGALModule from './sfcgal.js';

const Module = await SFCGALModule();
const sfcgal = new Module.SFCGAL();
sfcgal.initialize();

// Use sfcgal methods...

// Cleanup when done
sfcgal.delete();
```

### Available Operations

**Measurements:**
- `distance(wkt1, wkt2)` - Calculate 2D distance
- `distance3D(wkt1, wkt2)` - Calculate 3D distance
- `area(wkt)` - Calculate 2D area (validates geometry)
- `area3D(wkt)` - Calculate 3D area (validates geometry)
- `volume(wkt)` - Calculate volume (3D solids)
- `perimeter(wkt)` - Calculate perimeter (exterior ring for polygons)
- `perimeter3D(wkt)` - Calculate 3D perimeter

**2D Boolean Operations:**
- `intersection(wkt1, wkt2)` - 2D intersection
- `union(wkt1, wkt2)` - 2D union
- `difference(wkt1, wkt2)` - 2D difference

**3D Boolean Operations:**
- `intersection3D(wkt1, wkt2)` - 3D intersection (with validation and error handling)
- `union3D(wkt1, wkt2)` - 3D union (with validation and error handling)
- `difference3D(wkt1, wkt2)` - 3D difference (with validation and error handling)

**3D Operations:**
- `extrude(wkt, dx, dy, dz)` - Extrude geometry
- `buffer(wkt, distance)` - 2D buffer
- `buffer3D(wkt, distance)` - 3D buffer

**Validation:**
- `isValid(wkt)` - Check if geometry is valid
- `validity(wkt)` - Get detailed validity information

**Transformations:**
- `translate(wkt, dx, dy, dz)` - Translate geometry
- `scale(wkt, sx, sy, sz)` - Scale geometry
- `rotate(wkt, angle)` - Rotate in 2D
- `rotateX(wkt, angle)`, `rotateY(wkt, angle)`, `rotateZ(wkt, angle)` - 3D rotations

**Analysis:**
- `centroid(wkt)` - Calculate centroid
- `convexHull(wkt)` - Compute convex hull
- `orientation(wkt)` - Get orientation

**3D Primitives:**
- `makeSphere(x, y, z, radius, segments)` - Create sphere
- `makeCylinder(x, y, z, radius, height, segments)` - Create cylinder
- `makeCone(x, y, z, radius, height, segments)` - Create cone
- `makeTorus(x, y, z, majorRadius, minorRadius, majorSegments, minorSegments)` - Create torus

## Technical Details

### Dependencies

The binding compiles the following dependencies to WebAssembly:

- **GMP 6.3.0**: Arbitrary precision arithmetic
- **MPFR 4.2.2**: Multiple-precision floating-point
- **Boost 1.86.0**: C++ libraries (headers only)
- **CGAL 6.0.2**: Computational Geometry Algorithms Library
- **SFCGAL**: Latest version from `/home/lbartoletti/sfcgal`

### Build Configuration

- **Optimization**: `-O3` for maximum performance
- **Standards**: C++17
- **Memory**: 256MB initial, 2GB maximum, dynamic growth enabled
- **Exceptions**: Enabled (required for SFCGAL)
- **RTTI**: Enabled (required for Embind)
- **Special Flags**: `-DCGAL_DISABLE_ROUNDING_MATH_CHECK` (WebAssembly compatibility)

### File Sizes

- Uncompressed WASM: ~3.0 MB
- With gzip compression: ~800 KB (typical web server)
- With brotli compression: ~600 KB (optimal)

## Rendering Utilities

The examples include reusable rendering utilities:

### Canvas2D Renderer (`canvas2d-renderer.js`)

For 2D geometry visualization with HTML Canvas:

```javascript
import { Canvas2DRenderer } from './canvas2d-renderer.js';

const renderer = new Canvas2DRenderer('canvasId');

renderer.addGeometry('POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))', {
    strokeColor: '#667eea',
    fillColor: 'rgba(102, 126, 234, 0.2)',
    lineWidth: 2
});

renderer.render();
```

### Three3D Renderer (`three3d-renderer.js`)

For 3D geometry visualization with Three.js:

```javascript
import { Three3DRenderer } from './three3d-renderer.js';

const renderer = new Three3DRenderer('viewerId', {
    wireframe: false,
    showGrid: true,
    showAxes: true
});

renderer.addGeometry('POLYHEDRALSURFACE(...)', {
    color: 0x4caf50
});

// Toggle wireframe mode
renderer.setWireframe(true);
```

**Supported WKT Types:**
- POINT, LINESTRING, POLYGON, TRIANGLE
- POLYHEDRALSURFACE, TIN, SOLID
- GEOMETRYCOLLECTION (with nested geometries)

## License

This binding follows SFCGAL's license. SFCGAL is licensed under LGPL 2.0+.

## Credits

- **SFCGAL**: https://sfcgal.gitlab.io/SFCGAL/
- **CGAL**: https://www.cgal.org/
- **Emscripten**: https://emscripten.org/
- **Three.js**: Used in examples for 3D visualization
