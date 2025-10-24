file( GLOB SFCGAL_BENCHS_SOURCES *.cpp )
add_executable( bench-SFCGAL ${SFCGAL_BENCHS_SOURCES} )

find_package(Boost REQUIRED COMPONENTS unit_test_framework timer serialization)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

target_link_libraries( bench-SFCGAL
  SFCGAL
  ${Boost_LIBRARIES}
  ${GMP_LIBRARIES}
  ${MPFR_LIBRARIES}
)
set_target_properties( bench-SFCGAL PROPERTIES DEBUG_POSTFIX "d" )
install( TARGETS bench-SFCGAL DESTINATION bin )
