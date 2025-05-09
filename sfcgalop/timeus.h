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
 * @return Current time in microseconds, or -1.0 on error
 */
static inline double
get_time_us(void)
{
  static volatile LONG initialized = 0;
  static LARGE_INTEGER freq;

  /* Thread-safe initialization using InterlockedCompareExchange */
  if (InterlockedCompareExchange(&initialized, 1, 0) == 0) {
    if (!QueryPerformanceFrequency(&freq)) {
      InterlockedExchange(&initialized, 0);
      return -1.0;
    }

    /* Validate frequency to prevent division by zero */
    if (freq.QuadPart <= 0) {
      InterlockedExchange(&initialized, 0);
      return -1.0;
    }
  } else {
    /* Wait for initialization to complete */
    while (InterlockedCompareExchange(&initialized, 1, 1) == 0) {
      Sleep(0);
    }

    /* Check if initialization failed */
    if (freq.QuadPart <= 0) {
      return -1.0;
    }
  }

  LARGE_INTEGER counter;
  if (!QueryPerformanceCounter(&counter)) {
    return -1.0;
  }

  /* Explicit cast to prevent warnings and ensure correct conversion */
  const double counter_d = (double)(counter.QuadPart);
  const double freq_d    = (double)(freq.QuadPart);

  /* Check for potential overflow before multiplication */
  if (counter_d > (DBL_MAX / 1e6)) {
    return -1.0;
  }

  return counter_d * 1e6 / freq_d;
}

#elif defined(__APPLE__)
  #include <float.h>
  #include <mach/mach_time.h>

/**
 * Get current time in microseconds (Apple implementation)
 *
 * @return Current time in microseconds, or -1.0 on error
 */
static inline double
get_time_us(void)
{
  static mach_timebase_info_data_t timebase = {0, 0};
  static uint64_t                  start    = 0;

  /* Initialize timebase info once */
  if (timebase.numer == 0 || timebase.denom == 0) {
    kern_return_t kr = mach_timebase_info(&timebase);
    if (kr != KERN_SUCCESS) {
      return -1.0;
    }

    /* Validate timebase to prevent division by zero */
    if (timebase.denom == 0) {
      return -1.0;
    }

    start = mach_absolute_time();
    if (start == 0) {
      /* If mach_absolute_time returns 0, try again */
      start = mach_absolute_time();
    }
  }

  uint64_t now = mach_absolute_time();
  if (now < start) {
    /* Handle potential timer wrap-around */
    return -1.0;
  }

  uint64_t elapsed = now - start;

  /* Check for potential overflow in conversion */
  const double elapsed_d = (double)elapsed;
  const double numer_d   = (double)timebase.numer;
  const double denom_d   = (double)timebase.denom;

  if (elapsed_d > (DBL_MAX * denom_d / numer_d / 1000.0)) {
    return -1.0;
  }

  return elapsed_d * numer_d / denom_d / 1000.0;
}

#elif defined(__unix__) || defined(__linux__)
  #include <errno.h>
  #include <float.h>
  #include <time.h>

/**
 * Get current time in microseconds (Unix implementation)
 *
 * @return Current time in microseconds, or -1.0 on error
 */
static inline double
get_time_us(void)
{
  struct timespec ts;

  /* Use CLOCK_MONOTONIC for better reliability */
  if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) {
    /* Fallback to CLOCK_REALTIME if MONOTONIC fails */
    if (clock_gettime(CLOCK_REALTIME, &ts) != 0) {
      return -1.0;
    }
  }

  /* Validate timespec values */
  if (ts.tv_sec < 0 || ts.tv_nsec < 0 || ts.tv_nsec >= 1000000000L) {
    return -1.0;
  }

  /* Check for potential overflow */
  const double sec_us  = (double)ts.tv_sec * 1e6;
  const double nsec_us = (double)ts.tv_nsec / 1e3;

  if (sec_us > DBL_MAX - nsec_us) {
    return -1.0;
  }

  return sec_us + nsec_us;
}

#else
  #error                                                                       \
      "Unsupported platform for timeus.h - supported: Windows, Apple, Unix/Linux"
#endif

#endif /* TIMEUS_H */
