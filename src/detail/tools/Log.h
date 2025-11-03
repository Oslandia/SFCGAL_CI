// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_LOG_H_
#define SFCGAL_LOG_H_

#include "SFCGAL/config.h"

#include <boost/format.hpp>
#include <string>

/**
 *
 * Helper method to log debug message
 *
 * \code
 * SFCGAL_DEBUG( "start new method" ) ;
 * \endcode
 */
#define SFCGAL_DEBUG(message)                                                  \
  SFCGAL::Logger::get()->log(SFCGAL::Logger::Debug, message, __FILE__, __LINE__)
/**
 *
 * Helper method to log information message
 *
 * \code
 * SFCGAL_INFO( "start new method" ) ;
 * \endcode
 */
#define SFCGAL_INFO(message)                                                   \
  SFCGAL::Logger::get()->log(SFCGAL::Logger::Info, message, __FILE__, __LINE__)
/**
 *
 * Helper method to log warning message
 *
 * \code
 * SFCGAL_WARNING( "start new method" ) ;
 * \endcode
 */
#define SFCGAL_WARNING(message)                                                \
  SFCGAL::Logger::get()->log(SFCGAL::Logger::Warning, message, __FILE__,       \
                             __LINE__)
/**
 *
 * Helper method to log error message
 *
 * \code
 * SFCGAL_ERROR( "invalid geometry" ) ;
 * \endcode
 */
#define SFCGAL_ERROR(message)                                                  \
  SFCGAL::Logger::get()->log(SFCGAL::Logger::Info, message, __FILE__, __LINE__)
/**
 *
 * Helper method to log critical message
 *
 * \code
 * SFCGAL_ERROR( "unexpected behavior in triangulate" ) ;
 * \endcode
 */
#define SFCGAL_CRITICAL(message)                                               \
  SFCGAL::Logger::get()->log(SFCGAL::Logger::Critical, message, __FILE__,      \
                             __LINE__)

namespace SFCGAL {

/**
 * [Singleton]Logger class
 *
 * @warning saved_lines and co removed (dangerous for memory and could be done
 * in a LogWriter).
 */
class SFCGAL_API Logger {
public:
  /**
   * destructor
   */
  ~Logger();

  /**
   * log level
   */
  using Level = enum { Debug, Info, Warning, Error, Critical };

  /**
   * singleton accessor
   * @return Pointer to the global Logger instance
   */
  static Logger *
  get();

  /**
   * log a message using boost format
   * @param level the log level
   * @param message the message to log
   * @param filename the filename (optional)
   * @param lineNumber the line number in the file (optional)
   */
  void
  log(const Level &level, const boost::format &message,
      const std::string &filename = "", const int &lineNumber = -1);

  /**
   * log a message
   * @param level the log level
   * @param message the message to log
   * @param filename the filename (optional)
   * @param lineNumber the line number in the file (optional)
   */
  void
  log(const Level &level, const std::string &message,
      const std::string &filename = "", const int &lineNumber = -1);

  /**
   * get the current log level
   * @return The current logging level
   */
  const Level &
  logLevel() const;
  /**
   * set the log level
   * @param logLevel The new log level to set
   */
  void
  setLogLevel(const Level &logLevel);

private:
  /**
   * current log level
   */
  Level _logLevel = Warning;
  /**
   * display file position?
   */
  bool _displayFilePosition = true;

  /**
   * private constructor
   */
  Logger(std::ostream &);
  /**
   * no copy constructor
   */
  Logger(const Logger &other);

  std::ostream _out;
};

/**
 * @brief get the logger
 * @return Reference to the global logger instance
 */
SFCGAL_API auto
logger() -> Logger &;

} // namespace SFCGAL

#define SFCGAL_LOG(level, msg)                                                 \
  do {                                                                         \
    SFCGAL::Logger::get() << "[" << (level) << " " << __FILE__ << ":"          \
                          << __LINE__ << "] " << msg << std::endl;             \
  } while (0)

#ifndef NDEBUG
  #define LOG_DEBUG(msg)                                                       \
    do {                                                                       \
      SFCGAL_LOG("DEBUG", msg);                                                \
    } while (0)
#else
  #define LOG_DEBUG(msg)                                                       \
    do {                                                                       \
    } while (0)
#endif
#define LOG_NOTICE(msg)                                                        \
  do {                                                                         \
    SFCGAL_LOG("NOTICE", msg);                                                 \
  } while (0)
#define LOG_ERROR(msg)                                                         \
  do {                                                                         \
    SFCGAL_LOG("ERROR", msg);                                                  \
  } while (0)

#endif
