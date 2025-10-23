#!/bin/bash
set -e
set -o pipefail

echo "=== Building SFCGAL WebAssembly Module ==="
echo

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="${PROJECT_DIR}/build"
SRC_DIR="${PROJECT_DIR}/src"
EXAMPLES_DIR="${PROJECT_DIR}/examples"

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# Setup Emscripten environment
echo -e "${YELLOW}Setting up Emscripten environment...${NC}"
EMSDK_DIR="${PROJECT_DIR}/../build/emsdk"
if [ ! -d "${EMSDK_DIR}" ]; then
    EMSDK_DIR="${PROJECT_DIR}/build/emsdk"
fi

if [ -d "${EMSDK_DIR}" ]; then
    echo -e "${YELLOW}Using Emscripten from ${EMSDK_DIR}...${NC}"
    source "${EMSDK_DIR}/emsdk_env.sh"
else
    echo -e "${RED}Error: Emscripten not found at ${EMSDK_DIR}${NC}"
    echo -e "${RED}Please install Emscripten first${NC}"
    exit 1
fi

# Create output directory
OUTPUT_DIR="${EXAMPLES_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Compile the WASM module
echo -e "${YELLOW}Compiling SFCGAL binding...${NC}"
START_TIME=$(date +%s)

# Check if dependencies exist, if not build them
DEPS_DIR="${PROJECT_DIR}/deps"
DEPS_INSTALL="${DEPS_DIR}/install"

if [ ! -f "${DEPS_INSTALL}/lib/libgmp.a" ]; then
    echo -e "${YELLOW}Dependencies not found, building them first...${NC}"
    "${SCRIPT_DIR}/build-deps.sh"
    if [ $? -ne 0 ]; then
        echo -e "${RED}Failed to build dependencies${NC}"
        exit 1
    fi
fi

# Check for SFCGAL static library
SFCGAL_WASM_BUILD="${PROJECT_DIR}/build/sfcgal-wasm"

if [ -f "${SFCGAL_WASM_BUILD}/src/libSFCGAL.a" ]; then
    echo -e "${GREEN}Using SFCGAL from build${NC}"
    SFCGAL_WASM_BUILD="${SFCGAL_WASM_BUILD}"
else
    echo -e "${YELLOW}SFCGAL static library not found${NC}"
    echo -e "${YELLOW}Please run:${NC}"
    echo -e "${YELLOW}  ./scripts/build-sfcgal-wasm.sh${NC}"
    exit 1
fi

# Find SFCGAL installation paths
SFCGAL_INCLUDE="${SFCGAL_WASM_BUILD}/include"
SFCGAL_LIB="${SFCGAL_WASM_BUILD}/src"

echo "SFCGAL include: ${SFCGAL_INCLUDE}"
echo "SFCGAL lib: ${SFCGAL_LIB}"

# Find WebAssembly dependencies
BOOST_VERSION="1_86_0"
GMP_INCLUDE="${DEPS_INSTALL}/include"
GMP_LIB="${DEPS_INSTALL}/lib"
MPFR_INCLUDE="${DEPS_INSTALL}/include"
MPFR_LIB="${DEPS_INSTALL}/lib"
BOOST_INCLUDE="${DEPS_DIR}/boost_${BOOST_VERSION}"
CGAL_INCLUDE="${DEPS_INSTALL}/include"

em++ "${SRC_DIR}/sfcgal-binding.cpp" \
    --bind \
    -I"${SFCGAL_INCLUDE}" \
    -I"${CGAL_INCLUDE}" \
    -I"${BOOST_INCLUDE}" \
    -I"${GMP_INCLUDE}" \
    -I"${MPFR_INCLUDE}" \
    -L"${SFCGAL_LIB}" \
    -L"${GMP_LIB}" \
    -L"${MPFR_LIB}" \
    -lSFCGAL \
    -lgmp \
    -lmpfr \
    -DCGAL_DISABLE_ROUNDING_MATH_CHECK \
    -s WASM=1 \
    -s MODULARIZE=1 \
    -s EXPORT_ES6=1 \
    -s EXPORT_NAME="SFCGALModule" \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s INITIAL_MEMORY=256MB \
    -s MAXIMUM_MEMORY=2GB \
    -s EXPORTED_RUNTIME_METHODS='["ccall","cwrap"]' \
    -s ENVIRONMENT='web,worker,node' \
    -s SINGLE_FILE=0 \
    -O3 \
    -std=c++17 \
    -frtti \
    -fexceptions \
    -o "${OUTPUT_DIR}/sfcgal.js"

if [ $? -eq 0 ]; then
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))

    echo -e "${GREEN}✅ Build successful! (${DURATION} seconds)${NC}"
    echo
    echo "Output files:"
    ls -lah "${OUTPUT_DIR}"/sfcgal.* | grep -E "(\.js|\.wasm)$"
    echo
    echo -e "To test the module, run: ${YELLOW}./scripts/serve.sh${NC}"
else
    echo -e "${RED}❌ Build failed!${NC}"
    exit 1
fi
