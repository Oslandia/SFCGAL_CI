// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Exception.h"
#include <utility>

namespace SFCGAL {

Exception::Exception() noexcept : _message("unknown exception") {}

Exception::Exception(std::string const &message) noexcept : _message(message) {}

Exception::~Exception() noexcept = default;

const char *
Exception::what() const noexcept
{
  return _message.c_str();
}

std::string
Exception::diagnostic() const noexcept
{
  return boost::diagnostic_information(*this);
}

// Definitions of constructors and destructors for derived classes

GeometryInvalidityException::GeometryInvalidityException(
    std::string const &message) noexcept
    : Exception(message)
{
}

GeometryInvalidityException::~GeometryInvalidityException() noexcept = default;

NotImplementedException::NotImplementedException(
    std::string const &message) noexcept
    : Exception(message)
{
}

NotImplementedException::~NotImplementedException() noexcept = default;

InappropriateGeometryException::InappropriateGeometryException(
    std::string const &message) noexcept
    : Exception(message)
{
}

InappropriateGeometryException::~InappropriateGeometryException() noexcept =
    default;

NonFiniteValueException::NonFiniteValueException(
    std::string const &message) noexcept
    : Exception(message)
{
}

NonFiniteValueException::~NonFiniteValueException() noexcept = default;

WktParseException::WktParseException(std::string const &message) noexcept
    : Exception(message)
{
}

WktParseException::~WktParseException() noexcept = default;

} // namespace SFCGAL
