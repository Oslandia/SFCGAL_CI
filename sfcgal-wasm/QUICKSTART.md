# Quick Start Guide - SFCGAL WebAssembly

Get SFCGAL WebAssembly running in 5 minutes!

## Prerequisites

You need:
- Linux or macOS
- Git, CMake, curl, tar, xz
- ~2GB free disk space
- 15-30 minutes for first build

## Installation

### 1. Install Emscripten

```bash
cd /path/to/sfcgal-wasm
./scripts/install-emscripten.sh
source build/emsdk/emsdk_env.sh
```

### 2. Build Everything

**Quick method (No-Patch, Recommended):**
```bash
./scripts/build-all-no-patch.sh   # All-in-one script, clean build
```

**Manual steps (No-Patch):**
```bash
./scripts/build-deps-no-patch.sh  # ~15 minutes (no CGAL patches)
./scripts/build-sfcgal-wasm-no-patch.sh  # ~8 minutes
./scripts/build.sh                # ~20 seconds (auto-detects library)
```

**Legacy method (Patched):**
```bash
./scripts/build-deps.sh           # ~15 minutes (auto-patches CGAL CMake files)
./scripts/build-sfcgal-wasm.sh    # ~8 minutes
./scripts/build.sh                # ~20 seconds (auto-detects library)
```

### 3. Test It

```bash
./scripts/serve.sh
```

Open http://localhost:8000 in your browser.

## Usage Example

```javascript
// Using the shared module loader (recommended)
import { loadSFCGAL } from './sfcgal-module.js';

const sfcgal = await loadSFCGAL();

// Calculate distance
const dist = sfcgal.distance("POINT(0 0)", "POINT(3 4)");
console.log(dist); // 5.0

// Union of polygons
const result = sfcgal.union(
    "POLYGON((0 0,10 0,10 10,0 10,0 0))",
    "POLYGON((5 5,15 5,15 15,5 15,5 5))"
);

// Extrude to 3D
const solid = sfcgal.extrude("POLYGON((0 0,1 0,1 1,0 1,0 0))", 0, 0, 2);

// Validate geometry
if (!sfcgal.isValid(result)) {
    console.error("Invalid geometry!");
}
```

## Available Operations

**Measurements:**
- `distance()`, `distance3D()` - 2D/3D distance
- `area()`, `area3D()` - 2D/3D area (with validation)
- `perimeter()`, `perimeter3D()` - Perimeter calculation
- `volume()` - Volume of 3D solids

**Boolean Operations:**
- 2D: `intersection()`, `union()`, `difference()`
- 3D: `intersection3D()`, `union3D()`, `difference3D()`

**3D Operations:**
- `extrude()` - Extrude to 3D
- `buffer()`, `buffer3D()` - Buffer operations
- `makeSphere()`, `makeCylinder()`, `makeCone()`, `makeTorus()` - 3D primitives

**Other:**
- `isValid()`, `validity()` - Validation
- `translate()`, `scale()`, `rotate()`, `rotateX/Y/Z()` - Transformations
- `centroid()`, `convexHull()` - Analysis

## Troubleshooting

**"Emscripten not found"**
```bash
source build/emsdk/emsdk_env.sh
```

**"Dependencies not built"**
```bash
./scripts/build-deps.sh
```

**"Boost not found" or CMake errors**
```bash
rm -rf deps/install/ build/sfcgal-wasm
./scripts/build-deps.sh  # Automatically patches CGAL CMake
./scripts/build-sfcgal-wasm.sh
```

**Module doesn't load in browser**
- Use HTTP server (not file://)
- Run `./scripts/serve.sh`

## Interactive Examples

The project includes comprehensive examples organized by dimension:

**2D Examples** (using Canvas 2D):
- Area, distance, centroid calculations
- Boolean operations (intersection, union, difference)
- Convex hull, buffer, transformations

**3D Examples** (using Three.js):
- 3D extrusion with wireframe/solid modes
- 3D boolean operations (intersection3D, union3D, difference3D)
- 3D transformations and primitives

**Features:**
- Real-time visualization
- WKT input/output
- Example presets for each operation
- Wireframe toggle for 3D
- "Show Only Result" option for boolean operations

See `examples/LIMITATIONS.md` for WebAssembly-specific constraints and best practices.

## Next Steps

- Read [README.md](README.md) for full documentation
- Check [BUILDING.md](BUILDING.md) for detailed build instructions
- See [SUMMARY.md](SUMMARY.md) for technical details
- Browse `examples/` for more demos
- Review `examples/LIMITATIONS.md` for important usage tips

## File Sizes

- **sfcgal.wasm**: 3.0 MB uncompressed (~800KB gzipped)
- **sfcgal.js**: 107 KB uncompressed (~30KB gzipped)

## Support

For issues:
1. Check [BUILDING.md](BUILDING.md) troubleshooting section
2. Review error messages carefully
3. Try a clean rebuild: `rm -rf build/ deps/`

## License

LGPL 2.0+ (follows SFCGAL license)
