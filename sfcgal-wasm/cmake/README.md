# CMake Toolchain for SFCGAL WebAssembly (No-Patch Approach)

## Overview

This directory contains a CMake toolchain file that enables building SFCGAL for WebAssembly **without patching CGAL CMake files**.

## The Problem

CGAL's CMake configuration calls `find_package(Boost)` internally, which fails with Emscripten when using Boost in header-only mode because:
1. Boost's header-only distribution doesn't provide CMake config files
2. Emscripten's environment doesn't have a system Boost installation

## Two Solutions

### Solution 1: Patch CGAL CMake Files (Current Default)

**Scripts:**
- `scripts/build-deps.sh` - Patches CGAL after installation
- `scripts/build-sfcgal-wasm.sh` - Uses patched CGAL

**How it works:**
1. Build CGAL normally
2. Patch `CGAL_SetupBoost.cmake` to skip `find_package(Boost)`
3. Patch `CGAL_SetupCGAL_CoreDependencies.cmake` similarly
4. Build SFCGAL with patched CGAL

**Advantages:**
- ✅ Proven to work reliably
- ✅ Simple build process
- ✅ Currently used and tested

**Disadvantages:**
- ❌ Modifies CGAL installation files
- ❌ Patches may break with CGAL updates
- ❌ Less portable

### Solution 2: Custom FindBoost.cmake ✅ Validated

**Scripts:**
- `scripts/build-deps-no-patch.sh` - Builds deps without patches
- `scripts/build-sfcgal-wasm-no-patch.sh` - Uses custom FindBoost
- `cmake/FindBoost.cmake` - Custom Boost finder for header-only mode

**How it works:**
1. Build CGAL normally (no patches)
2. Add `cmake/` directory to CMAKE_MODULE_PATH
3. When CGAL calls `find_package(Boost)`, CMake finds our custom `FindBoost.cmake` first
4. Our FindBoost.cmake provides all needed Boost variables for header-only mode

**Advantages:**
- ✅ No modifications to CGAL installation
- ✅ Cleaner, more maintainable approach
- ✅ Better for CGAL version updates
- ✅ More portable

**Disadvantages:**
- ⚠️ Newer approach (less battle-tested than patches)
- ⚠️ Requires understanding of CMake module system

## Using the No-Patch Approach (Recommended)

### 1. Build Dependencies (No Patches)

```bash
./scripts/build-deps-no-patch.sh
```

This builds GMP, MPFR, Boost, and CGAL **without applying any patches**.

### 2. Build SFCGAL with Custom FindBoost

```bash
./scripts/build-sfcgal-wasm-no-patch.sh
```

This uses `-DCMAKE_MODULE_PATH` to inject our custom `FindBoost.cmake`.

### 3. Success! ✅

The build should complete with:
```
✅ SFCGAL static library built successfully!
✅ No CGAL patches were needed!
Library: build/sfcgal-wasm-no-patch/src/libSFCGAL.a
```

### 4. If It Fails (Fallback)

If the no-patch approach doesn't work on your system:

```bash
# Clean up
rm -rf build/sfcgal-wasm-no-patch

# Fall back to the patched method
./scripts/build-deps.sh          # This will apply patches
./scripts/build-sfcgal-wasm.sh   # Use patched CGAL
```

## Custom FindBoost.cmake Details

The `cmake/FindBoost.cmake` file intercepts CMake's `find_package(Boost)` calls:

1. **CMake Module Path Priority**: When `-DCMAKE_MODULE_PATH=cmake/` is set, CMake searches this directory first for Find*.cmake files

2. **Provides Boost Configuration**:
   ```cmake
   set(Boost_FOUND TRUE)
   set(Boost_VERSION "108600")  # 1.86.0
   set(Boost_INCLUDE_DIRS "${BOOST_INCLUDE_DIR}")
   set(Boost_LIBRARIES "")  # Header-only, no libs needed
   ```

3. **Satisfies CGAL Requirements**:
   - Version check: Sets Boost_VERSION to meet CGAL's minimum (>= 1.72)
   - Components: Marks all requested components as found
   - Paths: Provides include directories

4. **No System Search**: Prevents CMake from searching system paths for Boost

## Testing Both Approaches

You can test both approaches side-by-side:

```bash
# Patched approach (default)
./scripts/build-deps.sh
./scripts/build-sfcgal-wasm.sh
# Output: build/sfcgal-wasm/src/libSFCGAL.a

# Toolchain approach (experimental)
./scripts/build-deps-no-patch.sh
./scripts/build-sfcgal-wasm-no-patch.sh
# Output: build/sfcgal-wasm-no-patch/src/libSFCGAL.a
```

Both should produce identical libraries.

## Recommendations

**Current Status (2025):** The **no-patch approach** is now validated and working! ✅

**Recommended:** Use `build-deps-no-patch.sh` + `build-sfcgal-wasm-no-patch.sh` for clean builds.

**Legacy Support:** The patched scripts remain available as fallback for compatibility.

## Next Steps

Now that the no-patch approach is validated:

1. ✅ **Keep both methods** - Patched as legacy, no-patch as recommended
2. 🔄 **Monitor stability** - Test with different CGAL/CMake versions
3. 📝 **Update main docs** - Reference no-patch method in README/QUICKSTART
4. 🚀 **Consider making it default** - After more testing

## Contributing

Help improve the no-patch approach:
- Test with different CGAL versions (5.x, 6.x)
- Test on different platforms (Linux, macOS, WSL)
- Report any issues or edge cases
- Suggest improvements to `FindBoost.cmake`
