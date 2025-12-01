# CMakeLists.txt for sfcgalop - a CLI utility for SFCGAL operations

include(CheckCXXCompilerFlag)

set(SFCGALOP_SOURCES
    main.cpp
    operations/operations.cpp
    constructors.cpp
    io.cpp
    util.cpp
    text_ui.cpp
    error_handler.cpp
)

if(WIN32)
    list(APPEND SFCGALOP_SOURCES getopt_win.cpp)
endif()

add_executable(sfcgalop ${SFCGALOP_SOURCES})

target_include_directories(sfcgalop PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${SFCGAL_INCLUDE_DIRS}
)

target_link_libraries(sfcgalop PRIVATE
    SFCGAL
    ${CGAL_3RD_PARTY_LIBRARIES}
    Boost::thread
    Boost::serialization
)

if(WIN32)
    target_link_libraries(sfcgalop PRIVATE kernel32)
elseif(UNIX)
    find_package(Threads REQUIRED)
    target_link_libraries(sfcgalop PRIVATE Threads::Threads m)
endif()

install(TARGETS sfcgalop
        RUNTIME DESTINATION bin
        COMPONENT cli)

if(SFCGAL_BUILD_TESTS)
    add_subdirectory(test)
endif()
