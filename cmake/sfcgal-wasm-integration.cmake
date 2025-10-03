# SFCGAL WebAssembly Integration
# This file should be included from the main CMakeLists.txt

# Add option for WebAssembly build
option(SFCGAL_BUILD_WASM "Build WebAssembly binding (requires Emscripten)" OFF)

if(SFCGAL_BUILD_WASM)
    message(STATUS "WebAssembly binding build enabled")

    # Check if we're already using Emscripten
    if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
        message(FATAL_ERROR "Cannot build WASM binding when already using Emscripten toolchain. Use the standalone build script instead.")
    endif()

    # Add custom target that runs the build script
    add_custom_target(sfcgal-wasm
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/sfcgal-wasm/build.sh
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/sfcgal-wasm
        COMMENT "Building SFCGAL WebAssembly module (this may take a while on first run)"
        VERBATIM
    )

    # Installation of WASM files
    install(
        DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/sfcgal-wasm/build/output/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/share/sfcgal/wasm
        OPTIONAL  # Don't fail if build hasn't been run yet
        PATTERN "*.js"
        PATTERN "*.wasm"
        PATTERN "*.html"
    )

    # Add to 'all' target if requested
    option(SFCGAL_WASM_BUILD_WITH_ALL "Include WASM build in 'all' target" OFF)
    if(SFCGAL_WASM_BUILD_WITH_ALL)
        add_custom_target(sfcgal-wasm-all ALL DEPENDS sfcgal-wasm)
    endif()

    # Print instructions
    message(STATUS "To build WASM module: make sfcgal-wasm")
    message(STATUS "WASM files will be installed to: ${CMAKE_INSTALL_PREFIX}/share/sfcgal/wasm")
endif()