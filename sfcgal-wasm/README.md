# SFCGAL-Wasm

**SFCGAL compiled to WebAssembly** — High-performance 2D and 3D geometric operations in the browser and Node.js.

[![Build Status](https://github.com/SFCGAL/sfcgal-wasm/workflows/Build/badge.svg)](https://github.com/SFCGAL/sfcgal-wasm/actions)
[![License: LGPL-2.0](https://img.shields.io/badge/License-LGPL%202.0-blue.svg)](LICENSE)

---

## Overview

SFCGAL-Wasm is a WebAssembly port of [SFCGAL](https://github.com/SFCGAL/SFCGAL), a C++ library providing advanced 2D and 3D spatial operations. This port enables powerful geometric computations directly in web browsers and Node.js environments without server-side processing.

### Key Features

- **🌐 Browser & Node.js Support**: Run geometric operations anywhere JavaScript runs
- **📐 Comprehensive Operations**: Area, length, distance, buffer, intersection, union, difference, and more
- **🔺 3D Geometry Support**: Extrusion, 3D boolean operations, volume calculations
- **⚡ High Performance**: Uses CGAL with exact arithmetic for robust results
- **💾 WASM64 Memory Model**: Supports large geometries up to 16GB
- **🎨 Interactive Examples**: 18 ready-to-use examples with visualization

### Supported Geometry Types

- **2D**: POINT, LINESTRING, POLYGON, MULTIPOINT, MULTILINESTRING, MULTIPOLYGON
- **3D**: POINT Z, LINESTRING Z, POLYGON Z, POLYHEDRALSURFACE, SOLID
- **All** geometries support WKT (Well-Known Text) format

---

## Quick Start

### Prerequisites

- **Browser**: Chrome 119+, Firefox 121+, Safari 17+, or Edge 119+ (WASM64 support required)
- **Node.js** (optional): v20+ with `--experimental-wasm-memory64` flag

### Installation

#### Option 1: Clone and Build

```bash
# Clone the repository
git clone https://github.com/SFCGAL/sfcgal-wasm.git
cd sfcgal-wasm

# Build (requires cmake, git, curl, make, python3)
./build-complete.sh

# Serve examples locally
./scripts/serve.sh
# Open http://localhost:8888
```

**Build time**: 20-30 minutes (downloads and compiles GMP, MPFR, Boost, CGAL, SFCGAL)
**Disk space**: ~2GB for dependencies and build artifacts

#### Option 2: Use Pre-built Files

Download the latest release from [GitHub Releases](https://github.com/SFCGAL/sfcgal-wasm/releases):
- `sfcgal.js` (143 lines)
- `sfcgal.wasm` (4.1 MB)

---

## Usage

### In Browser

```html
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>SFCGAL Example</title>
</head>
<body>
    <script type="module">
        import { loadSFCGAL } from './lib/sfcgal-module.js';

        (async () => {
            // Load SFCGAL WebAssembly module
            const sfcgal = await loadSFCGAL();
            console.log('SFCGAL version:', sfcgal.version());

            // Calculate area of a polygon
            const wkt = "POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))";
            const area = sfcgal.area(wkt);
            console.log('Area:', area); // 100.0

            // Buffer operation
            const buffered = sfcgal.buffer(wkt, 2.0);
            console.log('Buffered geometry:', buffered);

            // 3D extrusion
            const extruded = sfcgal.extrude(wkt, 5.0);
            console.log('Extruded solid:', extruded);
        })();
    </script>
</body>
</html>
```

### In Node.js

```javascript
// test.mjs
import SFCGALModule from './sfcgal.js';

const Module = await SFCGALModule();
const sfcgal = new Module.SFCGAL();
sfcgal.initialize();

// Calculate distance between two points
const dist = sfcgal.distance(
    "POINT(0 0)",
    "POINT(3 4)"
);
console.log('Distance:', dist); // 5.0

// 3D union
const result = sfcgal.union3D(
    "POLYGON Z ((0 0 0, 4 0 0, 4 3 0, 0 3 0, 0 0 0))",
    "POLYGON Z ((2 1 0, 6 1 0, 6 4 0, 2 4 0, 2 1 0))"
);
console.log('Union:', result);
```

Run with:
```bash
node --experimental-wasm-memory64 test.mjs
```

---

## API Reference

### Initialization

#### `loadSFCGAL(options?): Promise<SFCGALInstance>`

Load and initialize the SFCGAL WebAssembly module.

**Parameters:**
- `options.timeout` (number): Timeout in milliseconds (default: 30000)
- `options.wasmPath` (string): Custom path to WASM file (default: '../sfcgal.wasm')

**Returns:** Promise resolving to SFCGAL instance

**Throws:** `SFCGALError` if initialization fails

```javascript
const sfcgal = await loadSFCGAL({ timeout: 10000 });
```

### Core Operations

#### Metrics

- **`area(wkt: string): number`** — Calculate 2D area
- **`area3D(wkt: string): number`** — Calculate 3D surface area
- **`volume(wkt: string): number`** — Calculate volume of 3D solid
- **`length(wkt: string): number`** — Calculate 2D length
- **`length3D(wkt: string): number`** — Calculate 3D length
- **`perimeter(wkt: string): number`** — Calculate 2D perimeter
- **`perimeter3D(wkt: string): number`** — Calculate 3D perimeter
- **`distance(wkt1: string, wkt2: string): number`** — Calculate 2D distance
- **`distance3D(wkt1: string, wkt2: string): number`** — Calculate 3D distance

#### Geometric Operations

- **`centroid(wkt: string): string`** — Calculate centroid
- **`convexHull(wkt: string): string`** — Calculate convex hull
- **`buffer(wkt: string, distance: number): string`** — Buffer geometry

#### Boolean Operations (2D)

- **`intersection(wkt1: string, wkt2: string): string`** — Calculate intersection
- **`union(wkt1: string, wkt2: string): string`** — Calculate union
- **`difference(wkt1: string, wkt2: string): string`** — Calculate difference

#### Boolean Operations (3D)

- **`intersection3D(wkt1: string, wkt2: string): string`** — Calculate 3D intersection
- **`union3D(wkt1: string, wkt2: string): string`** — Calculate 3D union
- **`difference3D(wkt1: string, wkt2: string): string`** — Calculate 3D difference

⚠️ **Note**: For planar geometries, use 2D operations instead of 3D. See [LIMITATIONS.md](examples/LIMITATIONS.md) for details.

#### 3D Operations

- **`extrude(wkt: string, height: number): string`** — Extrude 2D geometry to 3D
- **`extrudeDetailed(wkt: string, height: number, options: object): object`** — Extrude with metrics
- **`toSolid(wkt: string): string`** — Convert POLYHEDRALSURFACE to SOLID

#### Transformations

- **`translate(wkt: string, dx: number, dy: number, dz: number): string`** — Translate geometry
- **`rotate(wkt: string, angle: number): string`** — Rotate geometry (radians)
- **`scale(wkt: string, sx: number, sy: number, sz: number): string`** — Scale geometry

#### Validation

- **`isValid(wkt: string): boolean`** — Check geometry validity
- **`getGeometryInfo(wkt: string): object`** — Get geometry metadata

---

## Examples

The project includes 18 interactive examples demonstrating all major SFCGAL operations:

### 2D Examples
- Area & Length calculation
- Buffer operations
- Centroid computation
- Convex hull
- Boolean operations (intersection, union, difference)
- Distance calculations
- Geometric transformations

### 3D Examples
- 3D area and volume
- 3D boolean operations
- 3D distance
- Extrusion from 2D to 3D
- 3D transformations

### Utilities
- Geometry validation
- Tessellation
- Interactive playground

**Access examples**: Run `./scripts/serve.sh` and open http://localhost:8888

---

## Performance Considerations

### Memory

- **Initial memory**: 512 MB
- **Maximum memory**: 16 GB (WASM64)
- **Typical usage**: < 100 MB for simple geometries

### Optimization Tips

1. **Use 2D operations for planar geometries** — Significantly faster than 3D equivalents
2. **Validate geometries** — Use `isValid()` before complex operations
3. **Limit precision** — WKT with 1-2 decimal places is sufficient for most use cases
4. **Simplify geometries** — Reduce vertex count before expensive operations
5. **Monitor memory** — Large SOLID operations can require 1-2 GB

See [LIMITATIONS.md](examples/LIMITATIONS.md) for known constraints and best practices.

---

## Building from Source

### Prerequisites

- **CMake** 3.20+
- **Git**
- **curl**
- **tar**
- **make**
- **Python 3**
- **20-30 minutes** build time
- **~2GB** disk space

### Build Process

```bash
# Clone repository
git clone https://github.com/SFCGAL/sfcgal-wasm.git
cd sfcgal-wasm

# Full build (installs Emscripten, builds dependencies, compiles SFCGAL, creates Wasm module)
./build-complete.sh

# Clean build
./build-complete.sh --clean
```

The build script will:
1. Install Emscripten SDK (3.1.50)
2. Download and compile GMP 6.3.0, MPFR 4.2.2, Boost 1.86.0, CGAL 6.0.2
3. Compile SFCGAL as WebAssembly static library
4. Build final WebAssembly binding with Emscripten

**Output files**:
- `examples/sfcgal.js` (143 lines)
- `examples/sfcgal.wasm` (4.1 MB)

### Individual Build Steps

```bash
# Install Emscripten only
./scripts/install-emscripten.sh

# Build dependencies only
./scripts/build-deps.sh

# Build SFCGAL library only
./scripts/build-sfcgal-wasm.sh

# Build final Wasm module only
./scripts/build.sh
```

---

## Architecture

### Project Structure

```
sfcgal-wasm/
├── src/
│   └── sfcgal-binding.cpp      # C++ WebAssembly binding (Emscripten embind)
├── examples/
│   ├── lib/
│   │   ├── sfcgal-module.js    # JavaScript wrapper for WASM loading
│   │   ├── canvas2d-renderer.js # 2D Canvas renderer for geometries
│   │   └── three3d-renderer.js  # Three.js 3D renderer
│   ├── 2d/                     # 2D operation examples
│   ├── 3d/                     # 3D operation examples
│   ├── misc/                   # Utility examples
│   ├── index.html              # Example gallery
│   └── playground.html         # Interactive playground
├── scripts/
│   ├── build-complete.sh       # Master build orchestration
│   ├── install-emscripten.sh   # Emscripten SDK installer
│   ├── build-deps.sh           # Dependency builder
│   ├── build-sfcgal-wasm.sh    # SFCGAL WebAssembly compiler
│   ├── build.sh                # Final Wasm binding compiler
│   └── serve.sh                # Local HTTP server
├── cmake/
│   ├── EmscriptenToolchain.cmake # Custom Emscripten toolchain
│   ├── FindBoost.cmake          # Boost finder for header-only mode
│   └── FindThreads.cmake        # Threading stub for WebAssembly
└── deps/                       # Downloaded dependencies (git-ignored)
```

### Technology Stack

- **Core**: SFCGAL (C++) with CGAL, GMP, MPFR, Boost
- **Compilation**: Emscripten 3.1.50 with WASM64, Closure Compiler, LTO
- **JavaScript**: ES6 modules with async/await
- **Rendering**: Canvas 2D API + Three.js (r170+)
- **Build**: Bash scripts + CMake + Make

---

## Browser Compatibility

| Browser | Minimum Version | WASM64 Support |
|---------|----------------|----------------|
| Chrome | 119+ | ✅ |
| Firefox | 121+ | ✅ |
| Safari | 17+ | ✅ |
| Edge | 119+ | ✅ |
| Node.js | 20+ | ✅ (with flag) |

**Required features**:
- WebAssembly with 64-bit memory (`MEMORY64`)
- ES6 modules
- Dynamic imports (for Three.js)

---

## Known Limitations

### Memory Constraints
- **Maximum memory**: 16 GB (theoretical), ~4 GB (practical browser limit)
- **Complex 3D operations**: May fail with "Cannot allocate memory" error
- **SOLID geometries**: Require significant memory (1-2 GB for moderate complexity)

### 3D Boolean Operations
- **Planar geometries**: Use 2D algorithms (`intersection()`, `union()`, `difference()`) instead of 3D versions
- **CGAL Nef polyhedra**: Can fail with "Undecidable conversion" errors for certain inputs
- **Recommendation**: Keep geometries simple (< 100 vertices) for WebAssembly demonstrations

See [examples/LIMITATIONS.md](examples/LIMITATIONS.md) for comprehensive details.

---

## Contributing

Contributions are welcome! Please follow these guidelines:

1. **Fork** the repository
2. **Create** a feature branch (`git checkout -b feature/amazing-feature`)
3. **Commit** your changes (`git commit -m 'Add amazing feature'`)
4. **Push** to the branch (`git push origin feature/amazing-feature`)
5. **Open** a Pull Request

### Development Setup

```bash
# Clone and build
git clone https://github.com/SFCGAL/sfcgal-wasm.git
cd sfcgal-wasm
./build-complete.sh

# Run examples
./scripts/serve.sh

# Run tests (coming soon)
npm test
```

### Code Standards

- **JavaScript**: ES6+ with JSDoc comments
- **C++**: C++17 with RAII patterns
- **Formatting**: Follow .editorconfig rules
- **Testing**: Write unit tests for new features

---

## License

This project is licensed under the **GNU Lesser General Public License v2.0 (LGPL-2.0)** — see the [LICENSE](LICENSE) file for details.

SFCGAL-Wasm links to:
- **SFCGAL**: LGPL-2.0
- **CGAL**: GPL-3.0 / LGPL-3.0 (depending on components)
- **GMP**: LGPL-3.0 / GPL-2.0
- **MPFR**: LGPL-3.0
- **Boost**: Boost Software License 1.0

---

## Acknowledgments

- **SFCGAL Team**: Original C++ library and algorithms
- **CGAL Project**: Computational geometry algorithms
- **Emscripten**: WebAssembly compilation toolchain
- **Three.js**: 3D visualization library

---

## Support

- **Issues**: [GitHub Issues](https://github.com/SFCGAL/sfcgal-wasm/issues)
- **Documentation**: [examples/LIMITATIONS.md](examples/LIMITATIONS.md)
- **SFCGAL**: [Official website](http://sfcgal.org/)

---

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for version history.

---

**Made with ❤️ by the SFCGAL community**
