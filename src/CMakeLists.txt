file( GLOB_RECURSE SFCGAL_SOURCES "*.cpp" )
file( GLOB_RECURSE SFCGAL_HEADERS "*.h" )

if (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC )
    list(REMOVE_ITEM SFCGAL_SOURCES ${CMAKE_SOURCE_DIR}/src/algorithm/alphaShapes.cpp)
    list(REMOVE_ITEM SFCGAL_HEADERS ${CMAKE_SOURCE_DIR}/src/algorithm/alphaShapes.h)
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

target_link_libraries( SFCGAL CGAL::CGAL CGAL::CGAL_Core)

if( CGAL_USE_GMPXX )
target_link_libraries( SFCGAL gmpxx )
endif()

target_link_libraries( SFCGAL ${Boost_LIBRARIES} )

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
