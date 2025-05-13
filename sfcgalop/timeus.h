/**
 * timeus.h - Cross-platform microsecond precision timer
 */

#ifndef TIMEUS_H
#define TIMEUS_H

#include <stdint.h>

#if defined(_WIN32)
    #include <windows.h>
    /**
     * Get current time in microseconds (Windows implementation)
     * 
     * @return Current time in microseconds
     */
    static inline double get_time_us(void) {
        static LARGE_INTEGER freq;
        static int initialized = 0;
        
        if (!initialized) {
            QueryPerformanceFrequency(&freq);
            initialized = 1;
        }
        
        LARGE_INTEGER counter;
        QueryPerformanceCounter(&counter);
        
        return (double)(counter.QuadPart) * 1e6 / freq.QuadPart;
    }

#elif defined(__APPLE__)
    #include <mach/mach_time.h>
    /**
     * Get current time in microseconds (Apple implementation)
     * 
     * @return Current time in microseconds
     */
    static inline double get_time_us(void) {
        static mach_timebase_info_data_t timebase;
        static uint64_t start = 0;
        
        if (start == 0) {
            mach_timebase_info(&timebase);
            start = mach_absolute_time();
        }
        
        uint64_t now = mach_absolute_time() - start;
        return (double)now * timebase.numer / timebase.denom / 1000.0;
    }

#elif defined(__unix__) || defined(__linux__)
    #include <time.h>
    /**
     * Get current time in microseconds (Unix implementation)
     * 
     * @return Current time in microseconds
     */
    static inline double get_time_us(void) {
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        return (double)ts.tv_sec * 1e6 + (double)ts.tv_nsec / 1e3;
    }

#else
    #error "Unsupported platform for timeus.h"
#endif

#endif /* TIMEUS_H */
