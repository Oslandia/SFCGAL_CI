#include <stdio.h>

#if defined(_WIN32)
    #include <windows.h>
    double get_time_us() {
        static LARGE_INTEGER freq;
        static BOOL initialized = FALSE;
        if (!initialized) {
            QueryPerformanceFrequency(&freq);
            initialized = TRUE;
        }
        LARGE_INTEGER counter;
        QueryPerformanceCounter(&counter);
        return (double)(counter.QuadPart) * 1e6 / freq.QuadPart;
    }

#elif defined(__APPLE__)
    #include <mach/mach_time.h>
    double get_time_us() {
        static mach_timebase_info_data_t timebase;
        static uint64_t start = 0;
        if (start == 0) {
            mach_timebase_info(&timebase);
            start = mach_absolute_time();
        }
        uint64_t now = mach_absolute_time() - start;
        return (double)now * timebase.numer / timebase.denom / 1000.0;
    }

#elif defined(__unix__)
    #include <time.h>
    double get_time_us() {
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        return (double)ts.tv_sec * 1e6 + (double)ts.tv_nsec / 1e3;
    }

#else
    #error "Unsupported platform"
#endif
