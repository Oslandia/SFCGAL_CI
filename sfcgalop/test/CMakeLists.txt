add_executable(test_sfcgalop
    test_minimal.cpp
    ../io.cpp
    ../util.cpp
    ../error_handler.cpp
    ../text_ui.cpp
    ../constructors.cpp
    ../operations/operations.cpp
)

target_include_directories(test_sfcgalop PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${SFCGAL_INCLUDE_DIRS}
)

target_link_libraries(test_sfcgalop PRIVATE
    SFCGAL
    Boost::unit_test_framework
    ${CGAL_3RD_PARTY_LIBRARIES}
)

add_test(NAME sfcgalop_unit_tests COMMAND test_sfcgalop)
