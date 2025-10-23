# Stub FindThreads for WebAssembly
# WebAssembly doesn't have native threads support in the traditional sense
# This stub prevents CGAL from failing when looking for threads

set(CMAKE_THREAD_LIBS_INIT "")
set(CMAKE_USE_WIN32_THREADS_INIT 0)
set(CMAKE_USE_PTHREADS_INIT 0)
set(CMAKE_HP_PTHREADS_INIT 0)
set(Threads_FOUND TRUE)

# Create an interface library for threads
if(NOT TARGET Threads::Threads)
    add_library(Threads::Threads INTERFACE IMPORTED)
    set_target_properties(Threads::Threads PROPERTIES
        INTERFACE_COMPILE_OPTIONS ""
    )
endif()

message(STATUS "WebAssembly: Using stub Threads support (no native threads)")
