// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/Exception.h>

namespace SFCGAL {

///
///
///
Exception::Exception() throw() : _message("unknown exception") {}

///
///
///
Exception::Exception(std::string const &message) throw() : _message(message) {}

///
///
///
Exception::~Exception() throw() {}

///
///
///
const char *
Exception::what() const throw()
{
  return _message.c_str();
}

///
///
///
std::string
Exception::diagnostic() const throw()
{
  return boost::diagnostic_information(*this);
}

} // namespace SFCGAL
