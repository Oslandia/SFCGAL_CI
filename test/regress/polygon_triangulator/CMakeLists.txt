#-- polygon_triangulator test

if( SFCGAL_USE_STATIC_LIBS )
  add_definitions( "-DSFCGAL_USE_STATIC_LIBS" )
endif()

file( GLOB_RECURSE SFCGAL_REGRESS_POLYGON_TRIANGULATOR_TEST_SOURCES *.cpp )

set( REGRESS_NAME test-regress-polygon_triangulator )
add_executable( ${REGRESS_NAME} ${SFCGAL_REGRESS_POLYGON_TRIANGULATOR_TEST_SOURCES} )

find_package(Boost REQUIRED COMPONENTS ${SFCGAL_Boost_COMPONENTS})

target_link_libraries( ${REGRESS_NAME}	SFCGAL)
target_link_libraries( ${REGRESS_NAME} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES})


set_target_properties( ${REGRESS_NAME} PROPERTIES DEBUG_POSTFIX "d" )
