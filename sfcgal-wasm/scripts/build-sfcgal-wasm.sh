#!/bin/bash
set -e
set -o pipefail

# Function to get number of CPU cores (cross-platform)
get_num_cores() {
    if command -v nproc >/dev/null 2>&1; then
        nproc
    else
        sysctl -n hw.ncpu 2>/dev/null || echo 4
    fi
}

echo "=== Building SFCGAL Static Library for WebAssembly (No Patch Method) ==="
echo

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DEPS_DIR="${PROJECT_DIR}/deps"
DEPS_INSTALL="${DEPS_DIR}/install"
SFCGAL_WASM_BUILD="${PROJECT_DIR}/build/sfcgal-wasm"
TOOLCHAIN_FILE="${PROJECT_DIR}/cmake/EmscriptenToolchain.cmake"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# Setup Emscripten
EMSDK_DIR="${PROJECT_DIR}/build/emsdk"
if [ -d "${EMSDK_DIR}" ]; then
    echo -e "${YELLOW}Using Emscripten from ${EMSDK_DIR}...${NC}"
    source "${EMSDK_DIR}/emsdk_env.sh"
else
    echo -e "${RED}Error: Emscripten not found${NC}"
    exit 1
fi

# Check dependencies
if [ ! -f "${DEPS_INSTALL}/lib/libgmp.a" ]; then
    echo -e "${RED}Dependencies not built. Run ./scripts/build-deps.sh first${NC}"
    exit 1
fi

BOOST_VERSION="1_86_0"
BOOST_INCLUDE="${DEPS_DIR}/boost_${BOOST_VERSION}"
# Auto-detect SFCGAL source directory (parent of sfcgal-wasm)
SFCGAL_SRC="$(dirname "$PROJECT_DIR")"

mkdir -p "${SFCGAL_WASM_BUILD}"
cd "${SFCGAL_WASM_BUILD}"

echo -e "${YELLOW}Configuring SFCGAL for WebAssembly...${NC}"
echo -e "${YELLOW}Using custom FindBoost.cmake from: ${PROJECT_DIR}/cmake${NC}"

# Use CMAKE_MODULE_PATH to inject our custom FindBoost.cmake
# This makes CMake use our FindBoost before searching elsewhere
# This approach avoids the need to patch CGAL CMake files
EMCC_CFLAGS="-sMEMORY64=1" emcmake cmake "${SFCGAL_SRC}" \
    -DCMAKE_MODULE_PATH="${PROJECT_DIR}/cmake" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_POLICY_DEFAULT_CMP0167=NEW \
    -DCMAKE_CXX_FLAGS="-DCGAL_DISABLE_ROUNDING_MATH_CHECK -sMEMORY64=1" \
    -DCMAKE_C_FLAGS="-sMEMORY64=1" \
    -DSFCGAL_USE_STATIC_LIBS=ON \
    -DBUILD_SHARED_LIBS=OFF \
    -DSFCGAL_BUILD_TESTS=OFF \
    -DSFCGAL_BUILD_EXAMPLES=OFF \
    -DSFCGAL_BUILD_CLI=OFF \
    -DGMP_INCLUDE_DIR="${DEPS_INSTALL}/include" \
    -DGMP_LIBRARIES="${DEPS_INSTALL}/lib/libgmp.a" \
    -DMPFR_INCLUDE_DIR="${DEPS_INSTALL}/include" \
    -DMPFR_LIBRARIES="${DEPS_INSTALL}/lib/libmpfr.a" \
    -DBOOST_INCLUDE_DIR="${BOOST_INCLUDE}" \
    -DCGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE=TRUE \
    -DCGAL_HEADER_ONLY=TRUE \
    -DCGAL_DIR="${DEPS_INSTALL}/lib/cmake/CGAL" \
    -DTHREADS_PREFER_PTHREAD_FLAG=OFF \
    -DCMAKE_USE_PTHREADS_INIT=OFF \
    -DCMAKE_THREAD_LIBS_INIT="" \
    -DCMAKE_HAVE_THREADS_LIBRARY=OFF

if [ $? -ne 0 ]; then
    echo -e "${RED}CMake configuration failed${NC}"
    echo -e "${YELLOW}Falling back to patched method...${NC}"
    echo -e "${YELLOW}The toolchain approach may not work with this CGAL version${NC}"
    exit 1
fi

echo -e "${YELLOW}Building SFCGAL static library...${NC}"
emmake make -j$(get_num_cores)

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✅ SFCGAL static library built successfully!${NC}"
    echo -e "${GREEN}✅ No CGAL patches were needed!${NC}"
    echo "Library: ${SFCGAL_WASM_BUILD}/src/libSFCGAL.a"
    echo
    echo "You can now use this library with:"
    echo "  SFCGAL_LIB=\"${SFCGAL_WASM_BUILD}/src\""
    echo "  SFCGAL_INCLUDE=\"${SFCGAL_WASM_BUILD}/include\""
else
    echo -e "${RED}❌ Build failed!${NC}"
    exit 1
fi
