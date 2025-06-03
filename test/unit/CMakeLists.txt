#-- build unit tests

find_package(Boost REQUIRED COMPONENTS ${SFCGAL_Boost_COMPONENTS})

if( SFCGAL_USE_STATIC_LIBS )
  add_definitions( "-DSFCGAL_USE_STATIC_LIBS" )
endif()

set(SFCGAL_CONFIG_EXECUTABLE "${CMAKE_BINARY_DIR}/sfcgal-config")
if(WIN32)
  set(SFCGAL_CONFIG_EXECUTABLE "${CMAKE_BINARY_DIR}/sfcgal-config.exe")
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/SFCGAL/testsfcgalconfig.cpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/SFCGAL/testsfcgalconfig.cpp"
  @ONLY
)

file(GLOB_RECURSE SFCGAL_UNIT_TEST_SOURCES
     "${CMAKE_CURRENT_BINARY_DIR}/testsfcgalconfig.cpp"
     "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

if (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC )
    list(REMOVE_ITEM SFCGAL_UNIT_TEST_SOURCES ${CMAKE_SOURCE_DIR}/test/unit/SFCGAL/algorithm/AlphaShapesTest.cpp)
endif()

add_executable( unit-test-SFCGAL ${SFCGAL_UNIT_TEST_SOURCES} )
target_link_libraries( unit-test-SFCGAL SFCGAL)
target_link_libraries(unit-test-SFCGAL ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES})

set_target_properties( unit-test-SFCGAL PROPERTIES DEBUG_POSTFIX "d" )

if (CMAKE_BUILD_TYPE STREQUAL "Debug" )
    add_test( unit-testd ${CMAKE_CURRENT_BINARY_DIR}/unit-test-SFCGALd --auto_start_dbg=y --log_level=all)
else()
    add_test( unit-test ${CMAKE_CURRENT_BINARY_DIR}/unit-test-SFCGAL --auto_start_dbg=y --log_level=all)
endif()
