#!/bin/bash

echo "=== Building SFCGAL WebAssembly Dependencies (No Patch Method) ==="
echo "This version does NOT patch CGAL CMake files"
echo "Instead, it relies on the CMake toolchain file to configure Boost"
echo

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DEPS_DIR="${PROJECT_DIR}/deps"
BUILD_DIR="${DEPS_DIR}/build"
INSTALL_DIR="${DEPS_DIR}/install"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# Setup Emscripten
EMSDK_DIR="${PROJECT_DIR}/build/emsdk"
if [ ! -d "${EMSDK_DIR}" ]; then
    EMSDK_DIR="${PROJECT_DIR}/../build/emsdk"
fi

if [ -d "${EMSDK_DIR}" ]; then
    echo -e "${YELLOW}Using Emscripten from ${EMSDK_DIR}...${NC}"
    source "${EMSDK_DIR}/emsdk_env.sh"
else
    echo -e "${RED}Error: Emscripten not found${NC}"
    echo "Run ./scripts/install-emscripten.sh first"
    exit 1
fi

# Create directories
mkdir -p "${BUILD_DIR}"
mkdir -p "${INSTALL_DIR}"
cd "${DEPS_DIR}"

# Versions
GMP_VERSION="6.3.0"
MPFR_VERSION="4.2.2"
BOOST_VERSION="1.86.0"
BOOST_VERSION_UNDERSCORE="1_86_0"
CGAL_VERSION="6.0.2"

echo
echo "Building dependencies without CGAL patches:"
echo "- GMP ${GMP_VERSION}"
echo "- MPFR ${MPFR_VERSION}"
echo "- Boost ${BOOST_VERSION} (header-only)"
echo "- CGAL ${CGAL_VERSION}"
echo

# === GMP ===
echo -e "${YELLOW}Building GMP ${GMP_VERSION}...${NC}"
if [ ! -f "${INSTALL_DIR}/lib/libgmp.a" ]; then
    if [ ! -f "gmp-${GMP_VERSION}.tar.xz" ]; then
        curl -LO "https://gmplib.org/download/gmp/gmp-${GMP_VERSION}.tar.xz"
    fi

    tar xf "gmp-${GMP_VERSION}.tar.xz"
    cd "gmp-${GMP_VERSION}"

    emconfigure ./configure \
        --prefix="${INSTALL_DIR}" \
        --enable-static \
        --disable-shared \
        --disable-assembly \
	--host=wasm32 \
        CFLAGS="-O3"

    emmake make -j$(nproc)
    emmake make install

    cd ..
    echo -e "${GREEN}✅ GMP built${NC}"
else
    echo -e "${GREEN}✅ GMP already built${NC}"
fi

# === MPFR ===
echo -e "${YELLOW}Building MPFR ${MPFR_VERSION}...${NC}"
if [ ! -f "${INSTALL_DIR}/lib/libmpfr.a" ]; then
    if [ ! -f "mpfr-${MPFR_VERSION}.tar.xz" ]; then
        curl -LO "https://www.mpfr.org/mpfr-current/mpfr-${MPFR_VERSION}.tar.xz"
    fi

    tar xf "mpfr-${MPFR_VERSION}.tar.xz"
    cd "mpfr-${MPFR_VERSION}"

    emconfigure ./configure \
        --prefix="${INSTALL_DIR}" \
        --enable-static \
        --disable-shared \
        --with-gmp="${INSTALL_DIR}" \
        CFLAGS="-O3"

    emmake make -j$(nproc)
    emmake make install

    cd ..
    echo -e "${GREEN}✅ MPFR built${NC}"
else
    echo -e "${GREEN}✅ MPFR already built${NC}"
fi

# === Boost (Header-only) ===
echo -e "${YELLOW}Downloading Boost ${BOOST_VERSION} (header-only)...${NC}"
if [ ! -d "boost_${BOOST_VERSION_UNDERSCORE}" ]; then
    if [ ! -f "boost_${BOOST_VERSION_UNDERSCORE}.tar.gz" ]; then
        curl -LO "https://archives.boost.io/release/${BOOST_VERSION}/source/boost_${BOOST_VERSION_UNDERSCORE}.tar.gz"
    fi

    tar xf "boost_${BOOST_VERSION_UNDERSCORE}.tar.gz"
    echo -e "${GREEN}✅ Boost extracted${NC}"
else
    echo -e "${GREEN}✅ Boost already extracted${NC}"
fi

# === CGAL ===
echo -e "${YELLOW}Building CGAL ${CGAL_VERSION} (header-only)...${NC}"
if [ ! -d "${INSTALL_DIR}/lib/cmake/CGAL" ]; then
    if [ ! -f "CGAL-${CGAL_VERSION}.tar.xz" ]; then
        curl -LO "https://github.com/CGAL/cgal/releases/download/v${CGAL_VERSION}/CGAL-${CGAL_VERSION}.tar.xz"
    fi

    tar xf "CGAL-${CGAL_VERSION}.tar.xz"
    cd "CGAL-${CGAL_VERSION}"
    mkdir -p build
    cd build

    emcmake cmake .. \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
        -DCMAKE_POLICY_DEFAULT_CMP0167=NEW \
        -DCGAL_HEADER_ONLY=TRUE \
        -DWITH_CGAL_Qt5=OFF \
        -DWITH_CGAL_Qt6=OFF \
        -DWITH_CGAL_ImageIO=OFF \
        -DWITH_CGAL_Core=ON \
        -DBoost_INCLUDE_DIR="${DEPS_DIR}/boost_${BOOST_VERSION_UNDERSCORE}" \
        -DBoost_NO_BOOST_CMAKE=TRUE \
        -DBoost_NO_SYSTEM_PATHS=TRUE \
        -DGMP_INCLUDE_DIR="${INSTALL_DIR}/include" \
        -DGMP_LIBRARIES="${INSTALL_DIR}/lib/libgmp.a" \
        -DMPFR_INCLUDE_DIR="${INSTALL_DIR}/include" \
        -DMPFR_LIBRARIES="${INSTALL_DIR}/lib/libmpfr.a"

    emmake make install

    cd ../..
    echo -e "${GREEN}✅ CGAL installed${NC}"
else
    echo -e "${GREEN}✅ CGAL already installed${NC}"
fi

echo
echo -e "${GREEN}✅ All dependencies built!${NC}"
echo
echo "Dependencies installed in: ${INSTALL_DIR}"
echo "GMP: ${INSTALL_DIR}/lib/libgmp.a"
echo "MPFR: ${INSTALL_DIR}/lib/libmpfr.a"
echo "Boost: ${DEPS_DIR}/boost_${BOOST_VERSION_UNDERSCORE}"
echo "CGAL: ${INSTALL_DIR}/lib/cmake/CGAL"
echo
echo -e "${YELLOW}NOTE: CGAL CMake files were NOT patched${NC}"
echo -e "${YELLOW}The toolchain file (cmake/EmscriptenToolchain.cmake) will handle Boost configuration${NC}"
echo
echo "Next steps:"
echo "  ./scripts/build-sfcgal-wasm.sh"
echo "  ./scripts/build.sh"
