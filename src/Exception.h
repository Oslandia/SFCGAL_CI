// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_EXCEPTION_H_
#define SFCGAL_EXCEPTION_H_

#include "SFCGAL/config.h"
#include <boost/exception/all.hpp>
#include <boost/format.hpp>
#include <boost/throw_exception.hpp>
#include <string>

namespace SFCGAL {

/**
 * Base SFCGAL Exception
 *
 * \code
 * BOOST_THROW_EXCEPTION( Exception("invalid geometry") );
 * \endcode
 */
class SFCGAL_API Exception : public virtual boost::exception,
                             public virtual std::exception {
public:
  /**
   * @brief Default constructor
   */
  Exception() noexcept;
  /**
   * @brief Constructor with message
   * @param message The exception message
   */
  explicit Exception(std::string const &message) noexcept;
  /**
   * @brief Copy constructor
   */
  Exception(const Exception &) noexcept = default;
  /**
   * @brief Copy assignment operator
   * @return Reference to this exception
   */
  auto
  operator=(const Exception &) noexcept -> Exception & = default;
  /**
   * @brief Move constructor
   */
  Exception(Exception &&) noexcept = default;
  /**
   * @brief Move assignment operator
   * @return Reference to this exception
   */
  auto
  operator=(Exception &&) noexcept -> Exception & = default;
  /**
   * @brief Destructor
   */
  ~Exception() noexcept override;

  /**
   * returns the exception message
   * @return The exception message
   */
  auto
  what() const noexcept -> const char * override;

  /**
   * returns diagnostic information (file, line, etc.)
   * @return Diagnostic information string
   */
  auto
  diagnostic() const noexcept -> std::string;

protected:
  /// @brief Exception message
  std::string _message;
};

/**
 * SFCGAL Exception thrown when invalid geometries are found before entering an
 * algo
 */
class SFCGAL_API GeometryInvalidityException : public Exception {
public:
  /**
   * @brief Constructor with message
   * @param message The exception message
   */
  explicit GeometryInvalidityException(std::string const &message) noexcept;
  /**
   * @brief Copy constructor
   */
  GeometryInvalidityException(const GeometryInvalidityException &) noexcept =
      default;
  /**
   * @brief Copy assignment operator
   * @return Reference to this exception
   */
  auto
  operator=(const GeometryInvalidityException &) noexcept
      -> GeometryInvalidityException & = default;
  /**
   * @brief Move constructor
   */
  GeometryInvalidityException(GeometryInvalidityException &&) noexcept =
      default;
  /**
   * @brief Move assignment operator
   * @return Reference to this exception
   */
  auto
  operator=(GeometryInvalidityException &&) noexcept
      -> GeometryInvalidityException & = default;
  /**
   * @brief Destructor
   */
  ~GeometryInvalidityException() noexcept override;
};

/**
 * SFCGAL Exception thrown when a function is not implemented
 */
class SFCGAL_API NotImplementedException : public Exception {
public:
  /**
   * @brief Constructor with message
   * @param message The exception message
   */
  explicit NotImplementedException(std::string const &message) noexcept;
  /**
   * @brief Copy constructor
   */
  NotImplementedException(const NotImplementedException &) noexcept = default;
  /**
   * @brief Copy assignment operator
   * @return Reference to this exception
   */
  auto
  operator=(const NotImplementedException &) noexcept
      -> NotImplementedException & = default;
  /**
   * @brief Move constructor
   */
  NotImplementedException(NotImplementedException &&) noexcept = default;
  /**
   * @brief Move assignment operator
   * @return Reference to this exception
   */
  auto
  operator=(NotImplementedException &&) noexcept
      -> NotImplementedException & = default;
  /**
   * @brief Destructor
   */
  ~NotImplementedException() noexcept override;
};

/**
 * SFCGAL Exception thrown when geometry is inappropriate for a function
 */
class SFCGAL_API InappropriateGeometryException : public Exception {
public:
  /**
   * @brief Constructor with message
   * @param message The exception message
   */
  explicit InappropriateGeometryException(std::string const &message) noexcept;
  /**
   * @brief Copy constructor
   */
  InappropriateGeometryException(
      const InappropriateGeometryException &) noexcept = default;
  /**
   * @brief Copy assignment operator
   * @return Reference to this exception
   */
  auto
  operator=(const InappropriateGeometryException &) noexcept
      -> InappropriateGeometryException & = default;
  /**
   * @brief Move constructor
   */
  InappropriateGeometryException(InappropriateGeometryException &&) noexcept =
      default;
  /**
   * @brief Move assignment operator
   * @return Reference to this exception
   */
  auto
  operator=(InappropriateGeometryException &&) noexcept
      -> InappropriateGeometryException & = default;
  /**
   * @brief Destructor
   */
  ~InappropriateGeometryException() noexcept override;
};

/**
 * SFCGAL Exception thrown when non finite value is found
 */
class SFCGAL_API NonFiniteValueException : public Exception {
public:
  /**
   * @brief Constructor with message
   * @param message The exception message
   */
  explicit NonFiniteValueException(std::string const &message) noexcept;
  /**
   * @brief Copy constructor
   */
  NonFiniteValueException(const NonFiniteValueException &) noexcept = default;
  /**
   * @brief Copy assignment operator
   * @return Reference to this exception
   */
  auto
  operator=(const NonFiniteValueException &) noexcept
      -> NonFiniteValueException & = default;
  /**
   * @brief Move constructor
   */
  NonFiniteValueException(NonFiniteValueException &&) noexcept = default;
  /**
   * @brief Move assignment operator
   * @return Reference to this exception
   */
  auto
  operator=(NonFiniteValueException &&) noexcept
      -> NonFiniteValueException & = default;
  ~NonFiniteValueException() noexcept override;
};

/**
 * SFCGAL Exception thrown when parsing WKT
 */
class SFCGAL_API WktParseException : public Exception {
public:
  /**
   * @brief Constructor with message
   * @param message The exception message
   */
  explicit WktParseException(std::string const &message) noexcept;
  /**
   * @brief Copy constructor
   */
  WktParseException(const WktParseException &) noexcept = default;
  /**
   * @brief Copy assignment operator
   * @return Reference to this exception
   */
  auto
  operator=(const WktParseException &) noexcept
      -> WktParseException & = default;
  /**
   * @brief Move constructor
   */
  WktParseException(WktParseException &&) noexcept = default;
  /**
   * @brief Move assignment operator
   * @return Reference to this exception
   */
  auto
  operator=(WktParseException &&) noexcept -> WktParseException & = default;
  /**
   * @brief Destructor
   */
  ~WktParseException() noexcept override;
};

} // namespace SFCGAL

#endif // SFCGAL_EXCEPTION_H_
