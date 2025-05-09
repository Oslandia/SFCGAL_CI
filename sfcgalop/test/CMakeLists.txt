# Unit tests for SFCGALOP

# Find Boost.Test
find_package(Boost REQUIRED COMPONENTS unit_test_framework)

# Create test executable
add_executable(test_sfcgalop
    test_minimal.cpp
    ../io.cpp
    ../util.cpp
    ../error_handler.cpp
    ../text_ui.cpp
    ../constructors.cpp
    ../operations/operations.cpp
)

# Include directories
target_include_directories(test_sfcgalop PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${SFCGAL_INCLUDE_DIRS}
)

# Link libraries
target_link_libraries(test_sfcgalop
    SFCGAL
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    ${CGAL_LIBRARIES}
    ${GMP_LIBRARIES}
    ${MPFR_LIBRARIES}
)

# Set C++ standard
set_target_properties(test_sfcgalop PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

# Add test
add_test(NAME sfcgalop_unit_tests COMMAND test_sfcgalop)

# Enable testing
enable_testing()