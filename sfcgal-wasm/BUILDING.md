# Building SFCGAL WebAssembly

This document provides detailed instructions for building the SFCGAL WebAssembly binding from source.

## Overview

The build process consists of four main steps:

1. **Install Emscripten SDK** - WebAssembly compiler toolchain
2. **Build Dependencies** - Compile GMP, MPFR, Boost, and CGAL for WebAssembly
3. **Build SFCGAL Library** - Compile SFCGAL as a static WebAssembly library
4. **Build the Binding** - Create the final WebAssembly module with Embind

## Prerequisites

### System Requirements

- **Operating System**: Linux or macOS
- **Disk Space**: ~2GB for dependencies and build artifacts
- **Memory**: 4GB RAM minimum, 8GB recommended
- **Time**: Initial build takes 20-30 minutes

### Required Tools

```bash
# Debian/Ubuntu
sudo apt-get install build-essential cmake git curl tar xz-utils

# Fedora/RHEL
sudo dnf install gcc gcc-c++ cmake git curl tar xz

# macOS
xcode-select --install
brew install cmake
```

## Step 1: Install Emscripten

Emscripten is the toolchain that compiles C++ to WebAssembly.

### Automatic Installation

```bash
cd /home/lbartoletti/sfcgal/production-ready
./scripts/install-emscripten.sh
```

This script will:
- Clone the Emscripten SDK repository
- Install the latest stable version
- Activate it for use

### Manual Installation

If you prefer to install Emscripten manually:

```bash
cd /home/lbartoletti/sfcgal/production-ready
mkdir -p build
cd build

# Clone emsdk
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk

# Install and activate
./emsdk install latest
./emsdk activate latest

# Activate in current shell
source ./emsdk_env.sh
```

### Verify Installation

```bash
source build/emsdk/emsdk_env.sh
emcc --version
```

You should see output like:
```
emcc (Emscripten gcc/clang-like replacement) 3.1.x
```

## Step 2: Build Dependencies

This step compiles all required dependencies for WebAssembly.

### Run the Build Script

```bash
./scripts/build-deps.sh
```

### What Gets Built

1. **GMP 6.3.0** (~5 minutes)
   - GNU Multiple Precision Arithmetic Library
   - Arbitrary precision integers
   - Configured with `--disable-assembly` for WebAssembly

2. **MPFR 4.2.2** (~3 minutes)
   - Multiple-precision floating-point library
   - Depends on GMP

3. **Boost 1.86.0** (~1 minute)
   - Downloads headers only (header-only libraries)
   - No compilation needed
   - ~400MB download

4. **CGAL 6.0.2** (~5 minutes)
   - Computational Geometry Algorithms Library
   - Configured in header-only mode
   - Installed cmake configuration files

### Build Output

All dependencies are installed in:
```
deps/
├── boost_1_86_0/           # Boost headers
├── gmp-6.3.0/              # GMP source
├── mpfr-4.2.2/             # MPFR source
├── CGAL-6.0.2/             # CGAL source
└── install/                # Compiled libraries
    ├── include/            # Headers (GMP, MPFR, CGAL)
    └── lib/                # Static libraries (.a files)
        ├── libgmp.a
        ├── libmpfr.a
        └── cmake/CGAL/     # CGAL CMake config
```

### Troubleshooting Dependencies

**Download failures:**
```bash
# Clear and retry
rm -rf deps/gmp-*.tar.xz deps/mpfr-*.tar.xz deps/boost*.tar.gz deps/cgal.tar.xz
./scripts/build-deps.sh
```

**Build failures:**
```bash
# Full clean rebuild
rm -rf deps/
./scripts/build-deps.sh
```

## Step 3: Build SFCGAL Library

This step compiles SFCGAL itself as a static WebAssembly library.

### Patch SFCGAL CMake Files

SFCGAL's build system needs modifications to work with WebAssembly:

```bash
./scripts/patch-sfcgal-cmake.sh
```

This script:
- Comments out `find_package(Boost)` calls
- Adds manual Boost configuration for WebAssembly
- Creates backups (.bak files) of original files

### Build SFCGAL

```bash
./scripts/build-sfcgal-wasm.sh
```

### Build Configuration

The build uses these CMake options:

```cmake
-DCMAKE_BUILD_TYPE=Release              # Optimized build
-DSFCGAL_USE_STATIC_LIBS=ON            # Build static library
-DBUILD_SHARED_LIBS=OFF                 # No shared libraries
-DSFCGAL_BUILD_TESTS=OFF                # Skip tests
-DSFCGAL_BUILD_EXAMPLES=OFF             # Skip examples
-DSFCGAL_BUILD_CLI=OFF                  # Skip CLI tool
-DBoost_INCLUDE_DIR=...                 # Our Boost headers
-DGMP_INCLUDE_DIR=...                   # Our GMP
-DMPFR_INCLUDE_DIR=...                  # Our MPFR
-DCGAL_DIR=...                          # Our CGAL config
```

### Build Output

```
build/sfcgal-wasm/
├── src/
│   └── libSFCGAL.a         # Static library (~50MB)
└── include/
    └── SFCGAL/             # Copied headers
```

### Compilation Details

- **Files compiled**: ~120 C++ source files
- **Time**: 5-10 minutes with parallel build
- **Compiler**: em++ (Emscripten C++ compiler)
- **Standard**: C++17
- **Optimization**: -O3

### Troubleshooting SFCGAL Build

**"Boost not found" errors:**
```bash
# Reapply patches
./scripts/patch-sfcgal-cmake.sh

# Clean and rebuild
rm -rf build/sfcgal-wasm
./scripts/build-sfcgal-wasm.sh
```

**"fatal error: boost/config.hpp":**
```bash
# Check Boost is downloaded
ls deps/boost_1_86_0/boost/config.hpp

# If missing, rebuild dependencies
rm -rf deps/boost_*
./scripts/build-deps.sh
```

**Linker errors:**
```bash
# Ensure all dependencies built correctly
ls -lh deps/install/lib/libgmp.a    # Should exist
ls -lh deps/install/lib/libmpfr.a   # Should exist

# Rebuild SFCGAL
rm -rf build/sfcgal-wasm
./scripts/build-sfcgal-wasm.sh
```

## Step 4: Build the Binding

Final step: compile the C++ binding and link everything together.

### Run the Build

```bash
./scripts/build.sh
```

### What Happens

1. **Checks for dependencies** - Builds them if missing
2. **Checks for SFCGAL library** - Builds it if missing
3. **Compiles the binding** - Uses Embind to create JavaScript bindings
4. **Links everything** - Combines all static libraries
5. **Generates output** - Creates sfcgal.js and sfcgal.wasm

### Compilation Command

The script runs em++ with:

```bash
em++ src/sfcgal-binding.cpp \
    --bind \                              # Use Embind
    -I build/sfcgal-wasm/include \        # SFCGAL headers
    -I deps/install/include \             # CGAL, GMP, MPFR
    -I deps/boost_1_86_0 \                # Boost headers
    -L build/sfcgal-wasm/src \            # SFCGAL library
    -L deps/install/lib \                 # GMP, MPFR
    -lSFCGAL -lgmp -lmpfr \               # Link libraries
    -s WASM=1 \                           # Generate WebAssembly
    -s MODULARIZE=1 \                     # ES6 module
    -s EXPORT_ES6=1 \                     # ES6 export
    -s ALLOW_MEMORY_GROWTH=1 \            # Dynamic memory
    -O3 \                                 # Maximum optimization
    -std=c++17 \                          # C++17 standard
    -frtti \                              # Enable RTTI (for Embind)
    -fexceptions \                        # Enable exceptions
    -o examples/sfcgal.js
```

### Build Output

```
examples/
├── sfcgal.js       # JavaScript loader and API (~107KB)
└── sfcgal.wasm     # WebAssembly binary (~3.0MB)
```

### File Size Analysis

**Uncompressed:**
- sfcgal.wasm: ~3.0 MB
- sfcgal.js: ~107 KB
- **Total**: ~3.1 MB

**With gzip (typical web server):**
- sfcgal.wasm: ~800 KB
- sfcgal.js: ~30 KB
- **Total**: ~830 KB

**With brotli (optimal):**
- sfcgal.wasm: ~600 KB
- sfcgal.js: ~25 KB
- **Total**: ~625 KB

### Troubleshooting Final Build

**"undefined reference" errors:**
```bash
# Check SFCGAL library exists
ls -lh build/sfcgal-wasm/src/libSFCGAL.a

# Rebuild SFCGAL
rm -rf build/sfcgal-wasm
./scripts/build-sfcgal-wasm.sh
./scripts/build.sh
```

**Memory errors:**
```bash
# Increase initial memory in scripts/build.sh
-s INITIAL_MEMORY=67108864  # 64MB instead of 16MB
```

**"Module not found" in browser:**
- Ensure you're serving via HTTP, not file:// protocol
- Run: `./scripts/serve.sh`

## Complete Build from Scratch

To do a complete clean build:

```bash
# Clean everything
rm -rf build/ deps/ examples/sfcgal.js examples/sfcgal.wasm

# Build step by step
./scripts/install-emscripten.sh
source build/emsdk/emsdk_env.sh
./scripts/build-deps.sh
./scripts/patch-sfcgal-cmake.sh
./scripts/build-sfcgal-wasm.sh
./scripts/build.sh
```

## Incremental Builds

### Update SFCGAL source only

```bash
# SFCGAL source is at /home/lbartoletti/sfcgal
# After modifying SFCGAL source:
rm -rf build/sfcgal-wasm
./scripts/build-sfcgal-wasm.sh
./scripts/build.sh
```

### Update binding only

```bash
# After modifying src/sfcgal-binding.cpp:
./scripts/build.sh
# No need to rebuild SFCGAL or dependencies
```

### Update dependencies

```bash
# To update to new CGAL/GMP/MPFR versions:
# 1. Edit version numbers in scripts/build-deps.sh
# 2. Clean and rebuild:
rm -rf deps/ build/sfcgal-wasm
./scripts/build-deps.sh
./scripts/build-sfcgal-wasm.sh
./scripts/build.sh
```

## Optimization Tips

### Faster Builds

1. **Use ccache** (already enabled if installed):
```bash
sudo apt-get install ccache  # Debian/Ubuntu
```

2. **Parallel builds** (already used via `-j$(nproc)`):
- Speeds up compilation significantly
- Uses all CPU cores

3. **Incremental builds**:
- Only rebuild what changed
- Keep build/ directory

### Smaller Output

1. **Enable link-time optimization**:
```bash
# Add to scripts/build.sh:
-flto
```

2. **Strip debug info** (already done with `-O3`)

3. **Use closure compiler**:
```bash
# Add to scripts/build.sh:
-s CLOSURE_COMPILER=1
```

## Build Environment Variables

You can customize the build with environment variables:

```bash
# Use different SFCGAL source
export SFCGAL_SRC=/path/to/sfcgal
./scripts/build-sfcgal-wasm.sh

# Change output directory
export OUTPUT_DIR=/path/to/output
./scripts/build.sh

# Use specific Emscripten version
source build/emsdk/emsdk_env.sh
./scripts/build.sh
```

## Continuous Integration

For CI/CD pipelines:

```bash
#!/bin/bash
set -e  # Exit on error

# Install dependencies
./scripts/install-emscripten.sh
source build/emsdk/emsdk_env.sh

# Build everything
./scripts/build-deps.sh
./scripts/patch-sfcgal-cmake.sh
./scripts/build-sfcgal-wasm.sh
./scripts/build.sh

# Verify output
test -f examples/sfcgal.js
test -f examples/sfcgal.wasm
echo "Build successful!"
```

## Advanced Configuration

### Custom Emscripten Flags

Edit `scripts/build.sh` and add flags to the `em++` command:

```bash
-s ASSERTIONS=1                    # Enable runtime assertions
-s SAFE_HEAP=1                     # Detect memory issues
-s STACK_OVERFLOW_CHECK=2          # Stack overflow detection
-s DEMANGLE_SUPPORT=1              # Better error messages
-g                                 # Debug symbols
```

### Custom SFCGAL Options

Edit `scripts/build-sfcgal-wasm.sh` and modify CMake options:

```bash
-DSFCGAL_WITH_OSG=ON              # Enable OpenSceneGraph support
-DSFCGAL_BUILD_EXAMPLES=ON         # Build examples
-DCMAKE_BUILD_TYPE=Debug           # Debug build
```

## Troubleshooting Guide

### Common Issues

1. **Emscripten not found**
   - Solution: `source build/emsdk/emsdk_env.sh`

2. **Dependencies missing**
   - Solution: `./scripts/build-deps.sh`

3. **CMake errors about Boost**
   - Solution: `./scripts/patch-sfcgal-cmake.sh`

4. **Linking errors**
   - Solution: Rebuild SFCGAL library

5. **Module loading errors in browser**
   - Solution: Use HTTP server, not file:// protocol

### Getting Help

If you encounter issues:

1. Check error messages carefully
2. Review this document
3. Try a clean rebuild
4. Check SFCGAL and Emscripten documentation

## Build Statistics

Typical build times on modern hardware (8-core CPU, SSD):

- **Emscripten install**: 5 minutes
- **Dependencies build**: 15 minutes
- **SFCGAL build**: 8 minutes
- **Binding build**: 20 seconds
- **Total (first time)**: ~30 minutes

Incremental build times:
- **Binding only**: 20 seconds
- **SFCGAL only**: 8 minutes
- **Dependencies only**: 15 minutes

Disk space usage:
- **Source files**: 500 MB
- **Build artifacts**: 1.5 GB
- **Final output**: 3.1 MB
- **Total**: ~2 GB
