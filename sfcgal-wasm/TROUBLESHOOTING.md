# Troubleshooting Guide - SFCGAL WebAssembly

## Common Build Issues

### 1. "archive member is neither Wasm object file nor LLVM bitcode"

**Symptom:**
```
wasm-ld: warning: libSFCGAL.a: archive member 'xxx.cpp.obj' is neither Wasm object file nor LLVM bitcode
wasm-ld: error: undefined symbol: SFCGAL::io::readWkt(...)
```

**Cause:** SFCGAL library was compiled with system compiler (g++) instead of Emscripten (em++).

**Solution:**

1. Check which dependencies you used:
   ```bash
   ls -la deps/install/lib/cmake/CGAL/CGAL_SetupBoost.cmake.bak
   ```

   - If file exists: You used patched deps, but tried no-patch SFCGAL build
   - If file doesn't exist: You used no-patch deps correctly

2. Clean and rebuild SFCGAL:
   ```bash
   rm -rf build/sfcgal-wasm-no-patch
   ./scripts/build-sfcgal-wasm-no-patch.sh
   ```

3. If still failing, ensure Emscripten is sourced:
   ```bash
   source build/emsdk/emsdk_env.sh
   which emcc  # Should point to build/emsdk/upstream/emscripten/emcc
   ./scripts/build-sfcgal-wasm-no-patch.sh
   ```

### 2. "CGAL appears to be patched"

**Symptom:**
```
Error: CGAL appears to be patched!
For no-patch build, you must use:
  rm -rf deps/install/lib/cmake/CGAL
  ./scripts/build-deps-no-patch.sh
```

**Cause:** Trying to use no-patch build scripts with patched CGAL.

**Solution:**

You have two options:

**Option A: Switch to no-patch (clean build):**
```bash
rm -rf deps/install/lib/cmake/CGAL
./scripts/build-deps-no-patch.sh
./scripts/build-sfcgal-wasm-no-patch.sh
./scripts/build.sh
```

**Option B: Use patched method:**
```bash
# Dependencies already patched, just build SFCGAL
./scripts/build-sfcgal-wasm.sh
./scripts/build.sh
```

### 3. "Boost not found" with No-Patch Method

**Symptom:**
```
CMake Error: Could not find a package configuration file provided by "Boost"
```

**Cause:** CMAKE_MODULE_PATH not set correctly, or FindBoost.cmake missing.

**Solution:**

1. Verify FindBoost.cmake exists:
   ```bash
   ls cmake/FindBoost.cmake
   ```

2. If missing, you may have an incomplete repository. Clone again or download FindBoost.cmake.

3. Check build script uses CMAKE_MODULE_PATH:
   ```bash
   grep CMAKE_MODULE_PATH scripts/build-sfcgal-wasm-no-patch.sh
   ```
   Should show: `-DCMAKE_MODULE_PATH="${PROJECT_DIR}/cmake"`

### 4. Mixed Patched and No-Patch Builds

**Symptom:** Confusion about which method is being used, or builds failing intermittently.

**Cause:** Both `build/sfcgal-wasm/` and `build/sfcgal-wasm-no-patch/` exist, causing conflicts.

**Solution:**

Clean everything and pick one method:

```bash
# Clean all builds
rm -rf build/sfcgal-wasm build/sfcgal-wasm-no-patch
rm -f examples/sfcgal.js examples/sfcgal.wasm

# Choose one method:

# No-Patch (Recommended):
./scripts/build-all-no-patch.sh

# OR Patched (Legacy):
./scripts/build-deps.sh
./scripts/build-sfcgal-wasm.sh
./scripts/build.sh
```

### 5. "Dependencies not built"

**Symptom:**
```
Dependencies not built. Run ./scripts/build-deps-no-patch.sh first
```

**Cause:** GMP/MPFR/Boost/CGAL not compiled for WebAssembly.

**Solution:**

```bash
./scripts/build-deps-no-patch.sh
```

Or for patched method:
```bash
./scripts/build-deps.sh
```

### 6. Emscripten Not Found

**Symptom:**
```
Error: Emscripten not found at build/emsdk
```

**Solution:**

1. Install Emscripten:
   ```bash
   ./scripts/install-emscripten.sh
   ```

2. Source environment:
   ```bash
   source build/emsdk/emsdk_env.sh
   ```

3. Verify:
   ```bash
   which emcc
   emcc --version
   ```

## Build Method Decision Tree

```
Do you have existing builds?
├─ No → Use no-patch method (recommended)
│   └─ ./scripts/build-all-no-patch.sh
│
└─ Yes → Check which method you used
    ├─ ls build/sfcgal-wasm-no-patch/ exists?
    │   └─ Continue with no-patch
    │       ./scripts/build.sh (auto-detects)
    │
    └─ ls build/sfcgal-wasm/ exists?
        └─ Continue with patched
            ./scripts/build.sh (auto-detects)
```

## Debugging Tips

### Check Emscripten Toolchain

```bash
source build/emsdk/emsdk_env.sh
echo "Emcc: $(which emcc)"
echo "Em++: $(which em++)"
emcc --version
```

### Verify SFCGAL Library is WebAssembly

```bash
file build/sfcgal-wasm-no-patch/src/libSFCGAL.a
# Should mention "current ar archive" with LLVM bitcode
```

### Check Build Detection

```bash
ls -l build/sfcgal-wasm-no-patch/src/libSFCGAL.a
ls -l build/sfcgal-wasm/src/libSFCGAL.a

# build.sh will use the first one it finds (no-patch has priority)
```

### Verbose CMake Output

Edit `build-sfcgal-wasm-no-patch.sh` and add:
```bash
-DCMAKE_VERBOSE_MAKEFILE=ON
```

## Clean Build from Scratch

If all else fails, clean everything:

```bash
# Nuclear option - removes ALL builds and dependencies
rm -rf build/ deps/ examples/sfcgal.js examples/sfcgal.wasm

# Rebuild from scratch
./scripts/install-emscripten.sh
source build/emsdk/emsdk_env.sh
./scripts/build-all-no-patch.sh
```

## Getting Help

If you still have issues:

1. Check `cmake/README.md` for build method details
2. Review build logs carefully - look for compiler errors before linker errors
3. Ensure you're using consistent build method (don't mix patch/no-patch)
4. Verify Emscripten is properly sourced before each build step
