// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Exception.h"
#include <utility>

namespace SFCGAL {

Exception::Exception() noexcept : _message("unknown exception") {}

Exception::Exception(std::string const &message) noexcept
    : _message(std::move(message))
{
}

Exception::~Exception() noexcept = default;

auto
Exception::what() const noexcept -> const char *
{
  return _message.c_str();
}

auto
Exception::diagnostic() const noexcept -> std::string
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
