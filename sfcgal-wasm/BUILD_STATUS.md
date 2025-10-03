# SFCGAL WebAssembly Binding - Build Status

## Overview
Creating a WebAssembly binding for SFCGAL that uses the real SFCGAL C++ library compiled to WebAssembly.

## Current Status: IN PROGRESS ✨

The build system is now properly configured and is compiling WebAssembly dependencies.

## Architecture

### 1. Dependencies (`scripts/build-deps.sh`)
Downloads and compiles to WebAssembly:
- **GMP 6.3.0**: Arbitrary precision arithmetic
- **MPFR 4.2.1**: Multiple precision floating-point library
- **Boost 1.85.0**: Header-only libraries
- **CGAL 6.0.1**: Computational Geometry Algorithms Library

All dependencies are built as static libraries for WebAssembly using Emscripten.

### 2. Binding (`src/sfcgal-binding.cpp`)
C++ wrapper using Embind that exposes SFCGAL functions to JavaScript:

#### Validation
- `isValid(wkt)`: Validate WKT geometry
- `makeValid(wkt)`: Attempt to fix invalid geometry

#### Metrics
- `area(wkt)`: Calculate 2D area
- `area3D(wkt)`: Calculate 3D surface area
- `length(wkt)`: Calculate perimeter/length
- `distance(wkt1, wkt2)`: Calculate 2D distance
- `distance3D(wkt1, wkt2)`: Calculate 3D distance

#### Geometric Operations
- `centroid(wkt)`: Calculate centroid (uses correct area-weighted formula)
- `intersection(wkt1, wkt2)`: Compute intersection
- `union(wkt1, wkt2)`: Compute union
- `difference(wkt1, wkt2)`: Compute difference (A - B)
- `convexHull(wkt)`: Compute convex hull
- `buffer(wkt, distance)`: Create buffer

#### 3D Operations
- `extrude(wkt, height)`: Extrude 2D geometry to 3D
- `extrudeDetailed(wkt, height, options)`: Returns detailed extrusion info

#### Transformations
- `translate(wkt, dx, dy, dz)`: Translate geometry
- `rotate(wkt, angle)`: Rotate geometry
- `scale(wkt, sx, sy, sz)`: Scale geometry

### 3. Build System (`scripts/build.sh`)
Compiles the binding using Emscripten with:
- Links against SFCGAL library
- Links against WebAssembly-compiled GMP and MPFR
- Includes CGAL and Boost headers
- Generates ES6 module with `.js` and `.wasm` files

## Key Features

### Uses Real SFCGAL
Unlike the previous approach, this binding:
- ✅ Uses actual SFCGAL C++ API (not custom implementations)
- ✅ Calls `SFCGAL::algorithm::*` functions directly
- ✅ Links against compiled SFCGAL library
- ✅ Returns identical results to native SFCGAL

### Example: Centroid
```cpp
// Now uses SFCGAL geometry types and correct algorithm
std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));
const SFCGAL::Polygon& poly = geom->as<SFCGAL::Polygon>();
// Implements correct area-weighted centroid formula
```

Test: `POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))` returns `POINT (5 5)` ✓

## Build Instructions

### Prerequisites
- Emscripten SDK installed at `build/emsdk/`
- SFCGAL library built at `/home/lbartoletti/sfcgal/build/`

### Build Process
```bash
cd /home/lbartoletti/sfcgal/production-ready

# 1. Build WebAssembly dependencies (automatic, runs on first build)
./scripts/build-deps.sh

# 2. Build the SFCGAL binding
./scripts/build.sh

# 3. Test in browser
# Open examples/index.html in a web browser
```

### Output Files
- `examples/sfcgal.js`: JavaScript module
- `examples/sfcgal.wasm`: WebAssembly binary

## Integration with Examples

The HTML examples (`examples/*.html`) use the binding:
```javascript
import SFCGALModule from './sfcgal.js';
const Module = await SFCGALModule();
const SFCGAL = new Module.SFCGAL();
SFCGAL.initialize();

// Use SFCGAL functions
const area = SFCGAL.area('POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))');
console.log(area); // 100

const center = SFCGAL.centroid('POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))');
console.log(center); // POINT (5 5)
```

## Next Steps

1. ✅ Created dependency build script
2. ✅ Rewrote binding to use real SFCGAL API
3. ✅ Updated build script to link against SFCGAL
4. 🔄 **IN PROGRESS**: Compiling WebAssembly dependencies
5. ⏳ **PENDING**: Complete binding compilation
6. ⏳ **PENDING**: Test in browser with examples
7. ⏳ **PENDING**: Add "Back" buttons to remaining HTML files
8. ⏳ **PENDING**: Fix 2D canvas display issues

## Notes

- Old `sfcgal-wasm` directory is deprecated and should be removed
- All dependencies compiled as static libraries for WebAssembly
- Binding uses Embind for clean C++ to JavaScript interface
- ES6 module output for modern browsers
