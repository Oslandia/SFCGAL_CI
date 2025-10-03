#-- capi test

if( SFCGAL_USE_STATIC_LIBS )
  add_definitions( "-DSFCGAL_USE_STATIC_LIBS" )
endif()

file( GLOB_RECURSE SFCGAL_REGRESS_POLYGON_TRIANGULATOR_TEST_SOURCES *.c )

set( REGRESS_NAME test-regress-capi )
add_executable( ${REGRESS_NAME} ${SFCGAL_REGRESS_POLYGON_TRIANGULATOR_TEST_SOURCES} )

find_package(Boost REQUIRED COMPONENTS ${SFCGAL_Boost_COMPONENTS})

target_link_libraries( ${REGRESS_NAME} PUBLIC SFCGAL)
target_link_libraries( ${REGRESS_NAME} PRIVATE nlohmann_json::nlohmann_json)


set_target_properties( ${REGRESS_NAME} PROPERTIES DEBUG_POSTFIX "d" )
