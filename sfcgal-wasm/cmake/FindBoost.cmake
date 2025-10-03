# Custom FindBoost.cmake for SFCGAL WebAssembly build
# This file intercepts find_package(Boost) calls from CGAL
# and provides a pre-configured Boost setup for header-only mode

# This is a workaround to avoid patching CGAL's CMake files
# When CGAL calls find_package(Boost), CMake will find this file first
# if CMAKE_MODULE_PATH points to this directory

message(STATUS "Using custom FindBoost.cmake for WebAssembly (header-only mode)")

# Mark Boost as found
set(Boost_FOUND TRUE)

# Set version to satisfy CGAL's requirements (needs >= 1.72)
set(Boost_VERSION "108600")
set(Boost_VERSION_STRING "1.86.0")
set(Boost_MAJOR_VERSION "1")
set(Boost_MINOR_VERSION "86")
set(Boost_SUBMINOR_VERSION "0")

# Set include directories
if(DEFINED BOOST_INCLUDE_DIR)
    set(Boost_INCLUDE_DIRS "${BOOST_INCLUDE_DIR}")
    set(Boost_INCLUDE_DIR "${BOOST_INCLUDE_DIR}")
elseif(DEFINED Boost_INCLUDE_DIR)
    set(Boost_INCLUDE_DIRS "${Boost_INCLUDE_DIR}")
else()
    message(WARNING "BOOST_INCLUDE_DIR not set, Boost headers may not be found")
endif()

# No libraries needed for header-only mode
set(Boost_LIBRARIES "")
set(Boost_LIBRARY_DIRS "")

# Configuration flags
set(Boost_USE_STATIC_LIBS TRUE)
set(Boost_USE_MULTITHREADED TRUE)
set(Boost_NO_BOOST_CMAKE TRUE)
set(Boost_NO_SYSTEM_PATHS TRUE)

# Mark all components as found (header-only)
if(Boost_FIND_COMPONENTS)
    foreach(component ${Boost_FIND_COMPONENTS})
        string(TOUPPER ${component} UPPERCOMPONENT)
        set(Boost_${UPPERCOMPONENT}_FOUND TRUE)
        set(Boost_${UPPERCOMPONENT}_LIBRARY "")
    endforeach()
endif()

# Success message
message(STATUS "Boost ${Boost_VERSION_STRING} found (header-only)")
message(STATUS "Boost include dir: ${Boost_INCLUDE_DIRS}")

# Mark as advanced to clean up cmake-gui
mark_as_advanced(
    Boost_INCLUDE_DIR
    Boost_INCLUDE_DIRS
    Boost_LIBRARY_DIRS
)
