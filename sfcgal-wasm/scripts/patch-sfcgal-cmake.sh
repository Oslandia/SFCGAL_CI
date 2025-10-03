#!/bin/bash

echo "Patching SFCGAL CMakeLists for WebAssembly build..."

SFCGAL_SRC="/home/lbartoletti/sfcgal"

# Backup original files if not already backed up
if [ ! -f "${SFCGAL_SRC}/CMakeLists.txt.bak" ]; then
    cp "${SFCGAL_SRC}/CMakeLists.txt" "${SFCGAL_SRC}/CMakeLists.txt.bak"
fi

if [ ! -f "${SFCGAL_SRC}/sfcgalop/CMakeLists.txt.bak" ]; then
    cp "${SFCGAL_SRC}/sfcgalop/CMakeLists.txt" "${SFCGAL_SRC}/sfcgalop/CMakeLists.txt.bak"
fi

# Patch main CMakeLists.txt - comment out Boost find_package
sed -i '170s/^find_package( Boost/#find_package( Boost/' "${SFCGAL_SRC}/CMakeLists.txt"

# Add manual Boost setup after the commented line
sed -i '170a\
# Manual Boost setup for WebAssembly\
if(NOT Boost_FOUND)\
  set(Boost_FOUND TRUE)\
  set(Boost_INCLUDE_DIRS ${Boost_INCLUDE_DIR})\
  set(Boost_LIBRARIES "")\
  include_directories(SYSTEM ${Boost_INCLUDE_DIR})\
  message(STATUS "Using manual Boost configuration for WebAssembly")\
  message(STATUS "Boost_INCLUDE_DIR = ${Boost_INCLUDE_DIR}")\
endif()' "${SFCGAL_SRC}/CMakeLists.txt"

# Patch sfcgalop/CMakeLists.txt - comment out Boost find_package
sed -i '50s/^find_package(Boost/#find_package(Boost/' "${SFCGAL_SRC}/sfcgalop/CMakeLists.txt"

# Add manual Boost setup
sed -i '50a\
# Manual Boost setup for WebAssembly\
if(NOT Boost_FOUND)\
  set(Boost_FOUND TRUE)\
  message(STATUS "Using manual Boost configuration for WebAssembly (sfcgalop)")\
endif()' "${SFCGAL_SRC}/sfcgalop/CMakeLists.txt"

echo "✅ SFCGAL CMakeLists patched successfully"
echo "Backups saved as .bak files"
echo ""
echo "Note: CGAL patches are applied automatically by build-deps.sh"
