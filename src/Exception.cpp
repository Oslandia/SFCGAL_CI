// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/Exception.h>

#include <utility>

namespace SFCGAL {

///
///
///
Exception::Exception() noexcept : _message("unknown exception") {}

///
///
///
Exception::Exception(std::string message) noexcept : _message(std::move(message)) {}

///
///
///
Exception::~Exception() noexcept = default;

///
///
///
auto
Exception::what() const noexcept -> const char *
{
  return _message.c_str();
}

///
///
///
auto
Exception::diagnostic() const noexcept -> std::string
{
  return boost::diagnostic_information(*this);
}

} // namespace SFCGAL
