// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <iostream>

#include "SFCGAL/detail/tools/Log.h"

namespace SFCGAL {

Logger::~Logger() = default;

auto
Logger::get() -> Logger *
{
  static Logger log(std::cout);
  return &log;
}

void
Logger::log(const Level &level, const boost::format &message,
            const std::string &filename, const int &lineNumber)
{
  log(level, message.str(), filename, lineNumber);
}

void
Logger::log(const Level &level, const std::string &message,
            const std::string &filename, const int &lineNumber)
{
  if (level < _logLevel) {
    return;
  }

  // ptime now = second_clock::local_time();
  //_out << to_iso_string(now) << ":" ;

  if (_displayFilePosition && !filename.empty()) {
    _out << filename << ":";
  }

  if (_displayFilePosition && lineNumber >= 0) {
    _out << lineNumber << ":";
  }

  switch (level) {
  case Debug:
    _out << " debug: ";
    break;

  case Info:
    _out << " info: ";
    break;

  case Warning:
    _out << " warning: ";
    break;

  case Error:
    _out << " error: ";
    break;

  case Critical:
    _out << " critical: ";
    break;
  }

  _out << message << '\n';
}

auto
Logger::logLevel() const -> const Logger::Level &
{
  return _logLevel;
}

void
Logger::setLogLevel(const Level &logLevel)
{
  _logLevel = logLevel;
}

Logger::Logger(std::ostream &str) : _out(str.rdbuf()) {}

auto
logger() -> Logger &
{
  return *Logger::get();
}

} // namespace SFCGAL
