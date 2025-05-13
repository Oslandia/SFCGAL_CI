# CMakeLists.txt for sfcgalop - a CLI utility for SFCGAL operations

# Set minimum CMake version
cmake_minimum_required(VERSION 3.12)

# Project name
project(sfcgalop 
        VERSION 1.0.0
        DESCRIPTION "Command-line utility for SFCGAL operations"
        LANGUAGES C)

# Set C standard
set(CMAKE_C_STANDARD 23)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Sources for sfcgalop
set(SFCGALOP_SOURCES
    main.c
    operations/operations.c
    operations/operations_params.c
    operations/operations_registry.c
    operations/operations_validity.c
    operations/operations_metrics.c
    operations/operations_construction.c
    operations/operations_spatial_analysis.c
    operations/operations_transformation.c
    operations/operations_skeleton.c
    operations/operations_tesselation.c
    operations/operations_alphashapes.c
    operations/operations_visibility.c
    operations/operations_coordinate.c
    io.c
    util.c
)

# Target executable
add_executable(sfcgalop ${SFCGALOP_SOURCES})

# Include directories - ensure SFCGAL headers are available
target_include_directories(sfcgalop PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src/capi
    ${CMAKE_SOURCE_DIR}/include
)

# Link against SFCGAL library
target_link_libraries(sfcgalop PRIVATE SFCGAL)

# Compiler flags
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(sfcgalop PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
        -Wconversion
        -Wshadow
    )
elseif(MSVC)
    target_compile_options(sfcgalop PRIVATE
        /W4
        /WX
    )
endif()

# Installation
install(TARGETS sfcgalop
        RUNTIME DESTINATION bin
        COMPONENT cli)
