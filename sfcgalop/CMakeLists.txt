
# CMakeLists.txt for sfcgalop - a CLI utility for SFCGAL operations

# Set minimum CMake version
cmake_minimum_required(VERSION 3.12)

# Project name
project(sfcgalop
        DESCRIPTION "Command-line utility for SFCGAL operations"
        LANGUAGES CXX)

# Include CheckCXXCompilerFlag module for compiler flag checking
include(CheckCXXCompilerFlag)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
message(STATUS "Using C++17 standard")

#
# Cmake policies
#
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

set(SFCGALOP_SOURCES
    main.cpp
    operations/operations.cpp
    constructors.cpp
    io.cpp
    util.cpp
    text_ui.cpp
    error_handler.cpp
)

# Add Windows-specific getopt implementation
if(WIN32)
    list(APPEND SFCGALOP_SOURCES getopt_win.cpp)
endif()

# Target executable
add_executable(sfcgalop ${SFCGALOP_SOURCES})

# Include directories
target_include_directories(sfcgalop PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${SFCGAL_INCLUDE_DIRS}
)

# Find required dependencies
find_package(CGAL REQUIRED)
find_package(Boost REQUIRED COMPONENTS thread serialization)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

# Link against SFCGAL library and its dependencies
target_link_libraries(sfcgalop PRIVATE
    SFCGAL
    CGAL::CGAL
    ${Boost_LIBRARIES}
    ${GMP_LIBRARIES}
    ${MPFR_LIBRARIES}
)

# Platform-specific link libraries for thread-safe functions
if(WIN32)
    target_link_libraries(sfcgalop PRIVATE kernel32)
elseif(UNIX)
    find_package(Threads REQUIRED)
    target_link_libraries(sfcgalop PRIVATE Threads::Threads)

    target_link_libraries(sfcgalop PRIVATE m)
endif()

# Installation
install(TARGETS sfcgalop
        RUNTIME DESTINATION bin
        COMPONENT cli)

# Testing
if(SFCGAL_BUILD_TESTS)
    add_subdirectory(test)
endif()
