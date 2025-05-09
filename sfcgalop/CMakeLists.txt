# CMakeLists.txt for sfcgalop - a CLI utility for SFCGAL operations

# Set minimum CMake version
cmake_minimum_required(VERSION 3.12)

# Project name
project(sfcgalop 
        DESCRIPTION "Command-line utility for SFCGAL operations"
        LANGUAGES C)

# Set C standard
set(CMAKE_C_STANDARD 23)
set(CMAKE_C_STANDARD_REQUIRED ON)

option(ENABLE_SECURITY_HARDENING "Enable security hardening features" ON)
option(ENABLE_FORTIFY_SOURCE "Enable FORTIFY_SOURCE" ON)

set(SFCGALOP_SOURCES
    main.c
    safe_string.c
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

# Platform-specific link libraries for thread-safe functions
if(WIN32)
    target_link_libraries(sfcgalop PRIVATE kernel32)
elseif(UNIX)
    find_package(Threads REQUIRED)
    target_link_libraries(sfcgalop PRIVATE Threads::Threads)
    
    # Link librt for clock_gettime on older systems
    target_link_libraries(sfcgalop PRIVATE rt)
    
    target_link_libraries(sfcgalop PRIVATE m)
endif()

if(ENABLE_SECURITY_HARDENING)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(sfcgalop PRIVATE
            # Basic warnings
            -Wall
            -Wextra
            -Wpedantic
            -Werror
            -Wconversion
            -Wshadow
            
            # Security-specific warnings
            -Wformat-security
            -Wformat=2
            -Wundef
            -Wcast-align
            -Wcast-qual
            -Wwrite-strings
            -Wstrict-prototypes
            -Wmissing-prototypes
            -Wredundant-decls
            -Wnested-externs
            -Winline
            -Wno-long-long
            -Wunused-parameter
            -Wunused-variable
            -Wunused-function
            
            # Additional security flags
            -fstack-protector-strong
            -fstack-clash-protection
            -fcf-protection=full
        )
        
        target_link_options(sfcgalop PRIVATE
            -Wl,-z,relro
            -Wl,-z,now
            -Wl,-z,noexecstack
        )
        
        # Enable FORTIFY_SOURCE for release builds
        if(ENABLE_FORTIFY_SOURCE)
            target_compile_definitions(sfcgalop PRIVATE
                $<$<CONFIG:Release>:_FORTIFY_SOURCE=2>
                $<$<CONFIG:RelWithDebInfo>:_FORTIFY_SOURCE=2>
                $<$<CONFIG:MinSizeRel>:_FORTIFY_SOURCE=2>
            )
        endif()
        
    elseif(MSVC)
        target_compile_options(sfcgalop PRIVATE
            /W4                                 # High warning level
            /WX                                 # Treat warnings as errors
            /Wall                               # Enable all warnings
            /sdl                                # Enable additional security checks
            /guard:cf                           # Control Flow Guard
            /Qspectre                           # Spectre mitigation (if available)
        )
        
        # Security-related definitions for MSVC
        target_compile_definitions(sfcgalop PRIVATE
            _CRT_SECURE_NO_WARNINGS=0          # Don't disable security warnings
            _SCL_SECURE_NO_WARNINGS=0          # Don't disable STL security warnings
        )
        
        target_link_options(sfcgalop PRIVATE
            /DYNAMICBASE                        # ASLR support
            /NXCOMPAT                           # DEP support
            /GUARD:CF                           # Control Flow Guard
        )
    endif()
endif()

# Definitions for cross-platform compatibility
target_compile_definitions(sfcgalop PRIVATE
    # Enable POSIX functions on all platforms
    _POSIX_C_SOURCE=200809L
    
    # Platform-specific definitions
    $<$<PLATFORM_ID:Windows>:WIN32_LEAN_AND_MEAN>
    $<$<PLATFORM_ID:Windows>:NOMINMAX>
    
    # Security-related definitions
    SFCGALOP_SECURITY_HARDENED=1
)

# Debug vs Release specific settings
target_compile_definitions(sfcgalop PRIVATE
    $<$<CONFIG:Debug>:DEBUG=1>
    $<$<CONFIG:Debug>:SFCGALOP_DEBUG_MODE=1>
    $<$<NOT:$<CONFIG:Debug>>:NDEBUG=1>
)

# Optimization settings
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(sfcgalop PRIVATE -O2 -flto)
        target_link_options(sfcgalop PRIVATE -flto)
    elseif(MSVC)
        target_compile_options(sfcgalop PRIVATE /O2 /GL)
        target_link_options(sfcgalop PRIVATE /LTCG)
    endif()
endif()

# Installation
install(TARGETS sfcgalop
        RUNTIME DESTINATION bin
        COMPONENT cli)
