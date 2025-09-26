# CMakeLists.txt for sfcgalop - a CLI utility for SFCGAL operations

# Set minimum CMake version
cmake_minimum_required(VERSION 3.12)

# Project name
project(sfcgalop
        DESCRIPTION "Command-line utility for SFCGAL operations"
        LANGUAGES CXX)

# Include CheckCXXCompilerFlag module for compiler flag checking
include(CheckCXXCompilerFlag)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
message(STATUS "Using C++17 standard")

option(ENABLE_SECURITY_HARDENING "Enable security hardening features" ON)
option(ENABLE_FORTIFY_SOURCE "Enable FORTIFY_SOURCE" ON)

set(SFCGALOP_SOURCES
    main.cpp
    operations/operations.cpp
    constructors.cpp
    io.cpp
    util.cpp
    text_ui.cpp
    error_handler.cpp
)

# Add Windows-specific getopt implementation
if(WIN32)
    list(APPEND SFCGALOP_SOURCES getopt_win.cpp)
endif()

# Target executable
add_executable(sfcgalop ${SFCGALOP_SOURCES})

# Include directories - ensure SFCGAL headers are available
target_include_directories(sfcgalop PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src/capi
    ${CMAKE_SOURCE_DIR}/include
)

# Find required dependencies
find_package(CGAL REQUIRED)
find_package(Boost REQUIRED COMPONENTS thread serialization)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

# Link against SFCGAL library and its dependencies
target_link_libraries(sfcgalop PRIVATE
    SFCGAL
    CGAL::CGAL
    ${Boost_LIBRARIES}
    ${GMP_LIBRARIES}
    ${MPFR_LIBRARIES}
)

# Platform-specific link libraries for thread-safe functions
if(WIN32)
    target_link_libraries(sfcgalop PRIVATE kernel32)
elseif(UNIX)
    find_package(Threads REQUIRED)
    target_link_libraries(sfcgalop PRIVATE Threads::Threads)
    
    target_link_libraries(sfcgalop PRIVATE m)
endif()

if(ENABLE_SECURITY_HARDENING)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        # Check for availability of security flags
        check_cxx_compiler_flag("-fstack-clash-protection" HAVE_STACK_CLASH_PROTECTION)
        check_cxx_compiler_flag("-fcf-protection=full" HAVE_CF_PROTECTION)
        
        # Base warning and security flags
        set(SECURITY_COMPILE_FLAGS
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
            -Wredundant-decls
            -Wnested-externs
            -Winline
            -Wno-long-long
            -Wunused-parameter
            -Wunused-variable
            -Wunused-function
            
            # Always available security flags
            -fstack-protector-strong

            # Avoid CI errors on OpenBSD/mingw
            -Wno-unused-command-line-argument
        )
        
        # Add conditional security flags
        if(HAVE_STACK_CLASH_PROTECTION)
            # macOS may have stack clash protection, but not on arm64 aka m* chip
            EXECUTE_PROCESS( COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE )
            if ( NOT ${ARCHITECTURE} STREQUAL arm64 )
                list(APPEND SECURITY_COMPILE_FLAGS -fstack-clash-protection)
            endif()
        endif()
        
        if(HAVE_CF_PROTECTION)
            list(APPEND SECURITY_COMPILE_FLAGS -fcf-protection=full)
        endif()
        
        target_compile_options(sfcgalop PRIVATE ${SECURITY_COMPILE_FLAGS})
        
        if(CMAKE_SYSTEM_NAME MATCHES "Linux|FreeBSD|OpenBSD|NetBSD|DragonFly|SunOS")
            # ELF
            target_link_options(sfcgalop PRIVATE
                -Wl,-z,relro
                -Wl,-z,now
                -Wl,-z,noexecstack
            )
        elseif(APPLE)
            # troll^WmacOS Mach-O
            target_link_options(sfcgalop PRIVATE
                -Wl,-bind_at_load
                -Wl,-dead_strip
            )
        endif()

        
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
        )
        
        # Check for Spectre mitigation availability (newer MSVC versions)
        check_cxx_compiler_flag("/Qspectre" HAVE_SPECTRE_MITIGATION)
        if(HAVE_SPECTRE_MITIGATION)
            target_compile_options(sfcgalop PRIVATE /Qspectre)
        endif()
        
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
    _GNU_SOURCE
    _DEFAULT_SOURCE
    __XSI_VISIBLE
    
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
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(sfcgalop PRIVATE -O2)
        # Disabled LTO due to CGAL template instantiation issues
        # target_compile_options(sfcgalop PRIVATE -O2 -flto)
        # target_link_options(sfcgalop PRIVATE -flto)
    elseif(MSVC)
        target_compile_options(sfcgalop PRIVATE /O2 /GL)
        target_link_options(sfcgalop PRIVATE /LTCG)
    endif()
endif()

# Installation
install(TARGETS sfcgalop
        RUNTIME DESTINATION bin
        COMPONENT cli)

# Testing
option(BUILD_TESTING "Build unit tests" OFF)
if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()
