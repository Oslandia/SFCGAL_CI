#-- standalone test

if( SFCGAL_USE_STATIC_LIBS )
  add_definitions( "-DSFCGAL_USE_STATIC_LIBS" )
endif()

file( GLOB_RECURSE SFCGAL_REGRESS_STANDALONE_TEST_SOURCES *.cpp )
add_executable( standalone-regress-test-SFCGAL ${SFCGAL_REGRESS_STANDALONE_TEST_SOURCES} )

find_package(Boost REQUIRED COMPONENTS ${SFCGAL_Boost_COMPONENTS})

target_link_libraries( standalone-regress-test-SFCGAL SFCGAL ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES})

set_target_properties( standalone-regress-test-SFCGAL PROPERTIES DEBUG_POSTFIX "d" )

if (CMAKE_BUILD_TYPE STREQUAL "Debug" )
    add_test( standalone-regress-testd ${CMAKE_CURRENT_BINARY_DIR}/standalone-regress-test-SFCGALd --auto_start_dbg=y --log_level=all)
else()
    add_test( standalone-regress-test ${CMAKE_CURRENT_BINARY_DIR}/standalone-regress-test-SFCGAL --auto_start_dbg=y --log_level=all)
endif()
