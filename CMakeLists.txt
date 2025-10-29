cmake_minimum_required( VERSION 3.6 )

set( SFCGAL_VERSION_MAJOR 2 )
set( SFCGAL_VERSION_MINOR 3 )
set( SFCGAL_VERSION_PATCH 0 )

set( SFCGAL_VERSION "${SFCGAL_VERSION_MAJOR}.${SFCGAL_VERSION_MINOR}.${SFCGAL_VERSION_PATCH}" )
project( SFCGAL VERSION ${SFCGAL_VERSION} )

set( CMAKE_DEBUG_POSTFIX "d" )
set(CMAKE_CXX_STANDARD 17)

# ignore system installed SFCGAL to avoid conflicts during build/edit
list(APPEND CMAKE_IGNORE_PREFIX_PATH "/usr/include/SFCGAL")
list(APPEND CMAKE_IGNORE_PREFIX_PATH "/usr/local/include/SFCGAL")

#
# Cmake policies
#
if(POLICY CMP0077)
    cmake_policy(SET CMP0077 NEW)
endif()

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

#----------------------------------------------------------------------------
# build options
#----------------------------------------------------------------------------

# Windows M_PI
add_compile_definitions(_USE_MATH_DEFINES)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
endif()

option( SFCGAL_BUILD_EXAMPLES "build examples" OFF )
option( SFCGAL_BUILD_TESTS "build unit and regress tests" OFF )
option( SFCGAL_BUILD_BENCH "Build benchmarks" OFF )
option( SFCGAL_BUILD_DOC "Build documentation" OFF )

option( SFCGAL_WITH_OSG "Compile with OpenSceneGraph support" OFF )
option( SFCGAL_VALGRIND "Compile with Valgrind support" OFF )
option( SFCGAL_TSAN "Enable ThreadSanitizer (Debug only)" OFF)

# Configure CCache if available
if(NOT MSVC)
  option(USE_CCACHE "Use ccache" ON)
  if (USE_CCACHE)
    find_program(CCACHE_FOUND ccache)
    if(CCACHE_FOUND)
        set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
        set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
        message(STATUS "ccache found")
    endif(CCACHE_FOUND)
  endif(USE_CCACHE)
endif(NOT MSVC)

#-- include finders and co
set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules;${CMAKE_MODULE_PATH}" )

if (CMAKE_CXX_COMPILER MATCHES ".*clang")
    set(CMAKE_COMPILER_IS_CLANGXX 1)
endif ()

#-----------------------------------------------------------
# dependencies
#-----------------------------------------------------------

if( SFCGAL_VALGRIND )
  # valgrind tends to alter how CGAL works. To avoid the '-frounding-math' warning we need to activate this define:
  add_definitions( "-DCGAL_DISABLE_ROUNDING_MATH_CHECK" )
endif()

# Require Debug build if SFCGAL_TSAN is ON
if ( SFCGAL_TSAN )
    if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        message(FATAL_ERROR "SFCGAL_TSAN is enabled but CMAKE_BUILD_TYPE is not Debug. Please configure with -DCMAKE_BUILD_TYPE=Debug")
    endif()

    message(STATUS "ThreadSanitizer enabled (C & C++)")
    add_compile_options(-fsanitize=thread -g)
    add_link_options(-fsanitize=thread)
endif()

#-- find CGAL  ---------------------------------------------
if( CGAL_USE_GMPXX )
	add_definitions( "-DCGAL_USE_GMPXX" )
endif()
option( CGAL_USE_AUTOLINK "disable CGAL autolink" OFF )

if( ${CGAL_USE_AUTOLINK} )
	add_definitions( "-DCGAL_NO_AUTOLINK" )
endif()

find_package( CGAL 5.6 COMPONENTS Core REQUIRED )
message( STATUS "CGAL ${CGAL_VERSION} found" )

set(CGAL_3RD_PARTY_LIBRARIES CGAL::CGAL CGAL::CGAL_Core)

include_directories( ${CMAKE_BINARY_DIR}/include )

# --- copy public headers to build/include so clang-tidy finds them early ---
set(_sfcgal_hdr_src_dirs
  "${CMAKE_SOURCE_DIR}/include"
  "${CMAKE_SOURCE_DIR}/src"
)

foreach(_d IN LISTS _sfcgal_hdr_src_dirs)
  if(EXISTS "${_d}")
    file(GLOB_RECURSE _hdrs
      "${_d}/*.h"
      "${_d}/*.hpp"
      "${_d}/*.inl"
    )
    if(_hdrs)
      message(STATUS "Copying headers from ${_d} to ${CMAKE_BINARY_DIR}/include/SFCGAL")
      foreach(_h IN LISTS _hdrs)
	# preserve relative layout under build/include/SFCGAL
        file(RELATIVE_PATH _rel "${_d}" "${_h}")
	get_filename_component(_dest_dir "${CMAKE_BINARY_DIR}/include/SFCGAL/${_rel}" DIRECTORY)
        file(MAKE_DIRECTORY "${_dest_dir}")
        file(COPY "${_h}" DESTINATION "${_dest_dir}")
      endforeach()
    endif()
  endif()
endforeach()
# --------------------------------------------------------------------------


#-- BOOST --------------------------------------------------
option( Boost_USE_AUTO_LINK "boost use autolink" OFF )
set(Boost_USE_DEBUG_RUNTIME OFF)
# set(Boost_USE_STATIC_LIBS OFF)
option( Boost_USE_DEBUG_RUNTIME "Use Boost debug runtime" OFF )

if( NOT ${Boost_USE_AUTO_LINK} )
	add_definitions( "-DBOOST_ALL_NO_LIB" )
endif()

option( Boost_USE_STATIC_LIBS "boost use dynamic libraries" OFF )
if( Boost_USE_STATIC_LIBS )
	message( STATUS "Boost_USE_STATIC_LIBS=ON" )
	add_definitions( "-DBOOST_THREAD_USE_LIB" )
else()
	message( STATUS "Boost_USE_STATIC_LIBS=OFF" )
#	add_definitions( "-DBOOST_TEST_DYN_LINK" )
	add_definitions( "-DBOOST_ALL_DYN_LINK" )
endif()

option( Boost_USE_MULTITHREAD "boost use multithread libraries" ON )
if( ${Boost_USE_MULTITHREAD} )
	message( STATUS "Boost_USE_MULTITHREAD=ON" )
else()
	message( STATUS "Boost_USE_MULTITHREAD=OFF" )
endif()

#-- minimalist build allowed with boost version 1.69 or newer
set( SFCGAL_Boost_COMPONENTS thread serialization )
#-- program_options
if ( SFCGAL_BUILD_TESTS OR SFCGAL_BUILD_EXAMPLES OR SFCGAL_BUILD_OSG )
	set( SFCGAL_Boost_COMPONENTS unit_test_framework program_options ${SFCGAL_Boost_COMPONENTS} )
endif()
find_package( Boost 1.69.0 COMPONENTS ${SFCGAL_Boost_COMPONENTS} REQUIRED )

if( SFCGAL_USE_STATIC_LIBS )
  add_definitions( "-DSFCGAL_USE_STATIC_LIBS" )
  if (NOT MSVC)
    add_definitions( "-fPIC" )
  endif()
endif()

#-- nlohmann/json -----------------------------------------
find_package( nlohmann_json 3.11 REQUIRED )

#-- OpenScenegraph -----------------------------------------
if ( SFCGAL_WITH_OSG )
  find_package( OpenSceneGraph COMPONENTS osgDB osgUtil )
  if( ${OPENSCENEGRAPH_FOUND} )
    message( STATUS "OPENSCENEGRAPH_INCLUDE_DIRS = ${OPENSCENEGRAPH_INCLUDE_DIRS}" )
    message( STATUS "OPENSCENEGRAPH_LIBRARIES = ${OPENSCENEGRAPH_LIBRARIES}" )

    include_directories( SYSTEM ${OPENSCENEGRAPH_INCLUDE_DIRS} )
  endif()
endif()

#-- note that SYSTEM turns -I/path to -isystem and avoid warnings in CGAL and Boost
include_directories( SYSTEM
	${CGAL_INCLUDE_DIRS}
	${Boost_INCLUDE_DIRS}
)
link_directories(
	${CGAL_LIBRARY_DIRS}
	${Boost_LIBRARY_DIRS}
)

#-- Coverage ---------------------------------------------------------------------
option( SFCGAL_COVERAGE "Generate code coverage report" OFF)
if( SFCGAL_COVERAGE )
include(CodeCoverage)
append_coverage_compiler_flags()
set(COVERAGE_EXCLUDES
  '${PROJECT_SOURCE_DIR}/CGAL/*'
  '${PROJECT_SOURCE_DIR}/test/*')
setup_target_for_coverage_gcovr_xml(
     NAME sfcgal_coverage                    # New target name
     EXECUTABLE ctest -j ${PROCESSOR_COUNT} -R unit-test
     DEPENDENCIES SFCGAL unit-test-SFCGAL
     )
endif()

#-- Warnings, frounding-math and gprof  ------------------------------------------
option( SFCGAL_WARNING_AS_ERROR "fail the build on warnings" OFF )
option( SFCGAL_BUILD_WITH_GPROF "build with gprof" OFF )
if(MSVC)
	# Force to always compile with W4
	if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
		string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
	else()
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
	endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
	# Update if necessary
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-long-long -pedantic -Wpointer-arith -Wcast-align -Wcast-qual -Wno-overloaded-virtual -Wformat=2 -Winit-self -Wmissing-include-dirs -Wwrite-strings -Wno-error=undef")#-Wfloat-equal -Wconversion -Wshadow
    if( SFCGAL_WARNING_AS_ERROR )
		 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wno-error=format")
	endif()

	# Allows profiling with gprof
	if(SFCGAL_BUILD_WITH_GPROF)
		set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-pg")
	endif()
endif()


if(NOT DEFINED CMAKE_INSTALL_LIBDIR)
  set(_LIBDIR_DEFAULT "lib")
  # Override this default 'lib' with 'lib64' iff:
  #  - we are on Linux system but NOT cross-compiling
  #  - we are NOT on debian
  #  - we are on a 64 bits system
  # reason is: amd64 ABI: http://www.x86-64.org/documentation/abi.pdf
  # Note that the future of multi-arch handling may be even
  # more complicated than that: http://wiki.debian.org/Multiarch
  if(CMAKE_SYSTEM_NAME MATCHES "Linux"
      AND NOT CMAKE_CROSSCOMPILING
      AND NOT EXISTS "/etc/debian_version")
    if(NOT DEFINED CMAKE_SIZEOF_VOID_P)
      message(AUTHOR_WARNING
        "Unable to determine default CMAKE_INSTALL_LIBDIR directory because no target architecture is known. "
        "Please enable at least one language before including GNUInstallDirs.")
    else()
      if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
        set(_LIBDIR_DEFAULT "lib64")
      endif()
    endif()
  endif()
  set(CMAKE_INSTALL_LIBDIR "${_LIBDIR_DEFAULT}" CACHE PATH "object code libraries (${_LIBDIR_DEFAULT})")
endif()

if(NOT DEFINED CMAKE_INSTALL_FULL_LIBDIR)
    set(CMAKE_INSTALL_FULL_LIBDIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif()

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#SET(CMAKE_MACOSX_RPATH ON)

#-- generate library headers
configure_file( ${CMAKE_SOURCE_DIR}/src/config.h.cmake include/SFCGAL/config.h )
configure_file( ${CMAKE_SOURCE_DIR}/src/version.h.cmake include/SFCGAL/version.h )

enable_testing()

#note : not available on windows without export/import
OPTION( SFCGAL_USE_STATIC_LIBS "define if SFCGAL is build as a static library" OFF )

#-- build the library
add_subdirectory( src )

#-- build test (todo only if boost use dyn link)
add_subdirectory( test )

#-- build examples
if( SFCGAL_BUILD_EXAMPLES )
	add_subdirectory( example )
endif()

#-- doxygen documentation (allows make doc when doxygen is found)
if( SFCGAL_BUILD_DOC )
  add_subdirectory( doc )
endif()


#-- Option to build the command-line interface
option(SFCGAL_BUILD_CLI "Build the SFCGAL command-line interface tool" ON)
if(SFCGAL_BUILD_CLI)
    add_subdirectory(sfcgalop)
endif()

#-- install directories
install(DIRECTORY ${CMAKE_BINARY_DIR}/include DESTINATION .)


#-- sfcgal-config
if ( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )
  set( SFCGAL_LIB_NAME "SFCGAL${CMAKE_DEBUG_POSTFIX}" )
else()
  set( SFCGAL_LIB_NAME "SFCGAL" )
endif()

if ( UNIX )
  include(GNUInstallDirs)
endif()
#set( SFCGAL_LIB_NAME ${${CMAKE_BUILD_TYPE}
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/sfcgal-config.in ${CMAKE_CURRENT_BINARY_DIR}/sfcgal-config @ONLY)
install( PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/sfcgal-config DESTINATION bin )

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/sfcgal.pc.in ${CMAKE_CURRENT_BINARY_DIR}/sfcgal.pc @ONLY)
install( FILES ${CMAKE_CURRENT_BINARY_DIR}/sfcgal.pc DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(SFCGAL_INCLUDE_DIR "${CMAKE_INSTALL_INCLUDEDIR}")
set(SFCGAL_LIBRARY_DIR "${CMAKE_INSTALL_LIBDIR}")

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/SFCGALConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/SFCGALConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/SFCGAL"
    PATH_VARS SFCGAL_INCLUDE_DIR SFCGAL_LIBRARY_DIR
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/SFCGALConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/SFCGALConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/SFCGALConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/SFCGAL"
)

install(TARGETS SFCGAL
        EXPORT SFCGALTargets
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

install(EXPORT SFCGALTargets
        FILE SFCGALTargets.cmake
        NAMESPACE SFCGAL::
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/SFCGAL")
