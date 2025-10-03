# Emscripten Toolchain for SFCGAL WebAssembly
# This toolchain file configures CMake to build SFCGAL for WebAssembly
# without requiring patches to CGAL CMake files

# Platform identification
set(CMAKE_SYSTEM_NAME Emscripten)
set(CMAKE_SYSTEM_PROCESSOR x86)

# Pre-configure Boost to avoid CGAL's find_package(Boost) calls
# This is the key to avoiding CGAL CMake patches
set(Boost_FOUND TRUE CACHE BOOL "Boost found (header-only mode)")
set(Boost_INCLUDE_DIRS "${BOOST_INCLUDE_DIR}" CACHE PATH "Boost include directories")
set(Boost_LIBRARIES "" CACHE STRING "Boost libraries (none in header-only mode)")
set(Boost_NO_BOOST_CMAKE TRUE CACHE BOOL "Don't use Boost's CMake config")
set(Boost_NO_SYSTEM_PATHS TRUE CACHE BOOL "Don't search system paths for Boost")
set(Boost_USE_STATIC_LIBS TRUE CACHE BOOL "Use static Boost libraries")
set(Boost_USE_MULTITHREADED TRUE CACHE BOOL "Use multithreaded Boost")

# Mark Boost as already found to prevent re-searching
# This prevents CGAL from calling find_package(Boost) again
mark_as_advanced(FORCE Boost_DIR)

# Additional compiler flags for WebAssembly
set(CMAKE_CXX_FLAGS_INIT "-DCGAL_DISABLE_ROUNDING_MATH_CHECK")
set(CMAKE_C_FLAGS_INIT "-DCGAL_DISABLE_ROUNDING_MATH_CHECK")

# Enable verbose output for debugging if needed
# set(CMAKE_VERBOSE_MAKEFILE ON)
