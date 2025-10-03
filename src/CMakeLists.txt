file( GLOB_RECURSE SFCGAL_SOURCES "*.cpp" )
file( GLOB_RECURSE SFCGAL_HEADERS "*.h" )

if (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC )
    list(REMOVE_ITEM SFCGAL_SOURCES ${CMAKE_SOURCE_DIR}/src/algorithm/alphaShapes.cpp)
    list(REMOVE_ITEM SFCGAL_HEADERS ${CMAKE_SOURCE_DIR}/src/algorithm/alphaShapes.h)
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND
    CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 15.0 AND
    CGAL_VERSION VERSION_LESS 6.0.2)
  set(ALPHAWRAP_SRC ${CMAKE_SOURCE_DIR}/src/algorithm/alphaWrapping3D.cpp)
  list(FIND SFCGAL_SOURCES  ${ALPHAWRAP_SRC} ALPHA_WRAPPING_INDEX)
  if (ALPHA_WRAPPING_INDEX GREATER_EQUAL 0)
    # With inline enabled CGAL::alpha_wrap_3 may segfault with GCC 15
    # See: https://github.com/CGAL/cgal/issues/8876
    set_source_files_properties(${ALPHAWRAP_SRC} PROPERTIES COMPILE_OPTIONS "-fno-inline")
  endif()
endif()

file( GLOB_RECURSE SFCGAL_HEADERS_COPIED RELATIVE ${CMAKE_SOURCE_DIR}/src "*.h" )
add_custom_target(copy)
foreach (header ${SFCGAL_HEADERS_COPIED})
    add_custom_command(TARGET copy
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_SOURCE_DIR}/src/${header} ${CMAKE_BINARY_DIR}/include/SFCGAL/${header} DEPENDS ${CMAKE_SOURCE_DIR}/src/${header})
endforeach()

set( SFCGAL_OSG_HEADERS
  ${CMAKE_SOURCE_DIR}/src/io/osg.h
  ${CMAKE_SOURCE_DIR}/src/detail/io/OsgFactory.h )
set( SFCGAL_OSG_SOURCES
  ${CMAKE_SOURCE_DIR}/src/io/osg.cpp
  ${CMAKE_SOURCE_DIR}/src/detail/io/OsgFactory.cpp )
message( STATUS "removing OSG dependencies from the library")
list(REMOVE_ITEM SFCGAL_HEADERS ${SFCGAL_OSG_HEADERS})
list(REMOVE_ITEM SFCGAL_SOURCES ${SFCGAL_OSG_SOURCES})

if( SFCGAL_USE_STATIC_LIBS )
  add_definitions( "-DSFCGAL_USE_STATIC_LIBS" )
  if (NOT MSVC)
    add_definitions( "-fPIC" )
  endif()
  add_library(
    SFCGAL
    ${SFCGAL_HEADERS}
    ${SFCGAL_SOURCES}
    )

else()
  add_definitions( "-DSFCGAL_BUILD_SHARED" )
  add_library(
    SFCGAL SHARED
    ${SFCGAL_HEADERS}
    ${SFCGAL_SOURCES}
    )
endif()

add_dependencies( SFCGAL copy )

# set VERSION and SOVERSION
set_target_properties( SFCGAL PROPERTIES VERSION ${SFCGAL_VERSION}
                                         SOVERSION ${SFCGAL_VERSION_MAJOR} )

target_link_libraries( SFCGAL PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries( SFCGAL PRIVATE ${CGAL_3RD_PARTY_LIBRARIES})

if( CGAL_USE_GMPXX )
  target_link_libraries( SFCGAL PRIVATE gmpxx)
endif()

target_link_libraries( SFCGAL PRIVATE "${Boost_LIBRARIES}" )
set_target_properties(SFCGAL PROPERTIES
  INTERFACE_LINK_LIBRARIES ""
)

if ( ${Use_precompiled_headers} )
  if(PCHSupport_FOUND)
    # Add "-fPIC" for shared library build
    if( ${SFCGAL_USE_STATIC_LIBS} )
      set( pch_option "" )
    else()
      set( pch_option "-fPIC" )
    endif()
    add_precompiled_header(SFCGAL
      include/CGAL/all.h
      ${pch_option}
      )
  endif()
endif()

# install library
install(
	TARGETS
	SFCGAL
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	BUNDLE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# SFCGAL-osg
if ( SFCGAL_WITH_OSG )
  if( SFCGAL_USE_STATIC_LIBS )
    add_definitions( "-DSFCGAL_USE_STATIC_LIBS" )
    if (NOT MSVC)
      add_definitions( "-fPIC" )
    endif()
    add_library(
      SFCGAL-osg
      ${SFCGAL_OSG_HEADERS}
      ${SFCGAL_OSG_SOURCES}
      )

  else()
    add_definitions( "-DSFCGAL_BUILD_SHARED" )
    add_library(
      SFCGAL-osg SHARED
      ${SFCGAL_OSG_HEADERS}
      ${SFCGAL_OSG_SOURCES}
      )
  endif()
  target_link_libraries( SFCGAL-osg
    SFCGAL
    ${OPENSCENEGRAPH_LIBRARIES}
  )
  set_target_properties( SFCGAL-osg PROPERTIES VERSION ${SFCGAL_VERSION}
                                         SOVERSION ${SFCGAL_VERSION_MAJOR} )
  install(
	TARGETS
	SFCGAL-osg
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	BUNDLE DESTINATION ${CMAKE_INSTALL_LIBDIR}
   )
endif()
