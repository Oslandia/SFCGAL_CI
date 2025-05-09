/**
 * @file error_handler.hpp
 * @brief Error handling and validation utilities for SFCGALOP
 */

#ifndef SFCGALOP_ERROR_HANDLER_HPP
#define SFCGALOP_ERROR_HANDLER_HPP

#include <SFCGAL/Geometry.h>
#include <cstdint>
#include <exception>
#include <string>

/**
 * @namespace ErrorHandler
 * @brief Error handling and validation functionality
 */
namespace ErrorHandler {

/**
 * @enum ErrorCode
 * @brief Error codes for different failure scenarios
 */
enum class ErrorCode : std::uint8_t {
  SUCCESS               = 0, ///< Operation completed successfully
  INVALID_ARGUMENTS     = 1, ///< Invalid command line arguments
  INVALID_GEOMETRY      = 2, ///< Invalid geometry input
  OPERATION_FAILED      = 3, ///< Geometry operation failed
  IO_ERROR              = 4, ///< Input/output error
  MEMORY_ERROR          = 5, ///< Memory allocation failure
  UNSUPPORTED_OPERATION = 6  ///< Requested operation not supported
};

/**
 * @class SfcgalopException
 * @brief Custom exception class for SFCGALOP errors
 */
class SfcgalopException : public std::exception {
private:
  std::string message;
  ErrorCode   code;

public:
  /**
   * @brief Construct an SfcgalopException with a message and associated error
   * code
   * @param msg Human-readable error message describing the failure
   * @param error_code ErrorCode categorizing the error
   */
  SfcgalopException(std::string msg, ErrorCode error_code)
      : message(std::move(msg)), code(error_code)
  {
  }

  /**
   * @brief Returns the explanatory string for this exception
   * @return Null-terminated message describing the error
   * @note The returned pointer refers to the internal message storage and
   * remains valid only while the exception object exists and is not modified.
   */
  [[nodiscard]] const char *
  what() const noexcept override // NOLINT(modernize-use-trailing-return-type)
  {
    return message.c_str();
  }

  /**
   * @brief Returns the error code associated with this exception
   * @return The ErrorCode stored in the exception
   */
  [[nodiscard]] ErrorCode
  error_code() const // NOLINT(modernize-use-trailing-return-type)
  {
    return code;
  }
};

/**
 * @struct ValidationResult
 * @brief Result of geometry validation
 */
struct ValidationResult {
  bool        valid;   ///< Whether geometry is valid
  std::string reason;  ///< Reason for invalidity
  std::string details; ///< Additional validation details
};

/**
 * @brief Validate a single geometry
 * @param geom Geometry to validate
 * @return Validation result with details
 */
auto
validate_geometry(const SFCGAL::Geometry &geom) -> ValidationResult;

/**
 * @brief Validate inputs for a specific operation
 * @param operation Operation name
 * @param first_geom First geometry (required)
 * @param second_geom Second geometry (may be nullptr)
 * @return Validation result
 */
auto
validate_operation_inputs(const std::string      &operation,
                          const SFCGAL::Geometry *first_geom,
                          const SFCGAL::Geometry *second_geom)
    -> ValidationResult;

/**
 * @brief Report an error to stderr with formatting
 * @param context Context where error occurred
 * @param exception Exception that was thrown
 */
void
report_error(const std::string &context, const std::exception &exception);

/**
 * @brief Report validation error with details
 * @param result Validation result to report
 */
void
report_validation_error(const ValidationResult &result);

/**
 * @brief Print detailed geometry information
 * @param geom Geometry to describe
 */
void
print_geometry_info(const SFCGAL::Geometry &geom);

} // namespace ErrorHandler

#endif // SFCGALOP_ERROR_HANDLER_HPP