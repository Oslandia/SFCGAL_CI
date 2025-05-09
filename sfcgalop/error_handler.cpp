// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "error_handler.hpp"
#include "operations/operations.hpp"
#include <SFCGAL/Envelope.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/Validity.h>
#include <SFCGAL/algorithm/isSimple.h>
#include <SFCGAL/algorithm/isValid.h>
#include <iostream>
#include <sstream>

namespace ErrorHandler {

auto
validate_geometry(const SFCGAL::Geometry &geom) -> ValidationResult
{
  ValidationResult result{true, "", ""};

  // Check validity and get detailed error message
  SFCGAL::Validity validity = SFCGAL::algorithm::isValid(geom);

  if (!validity) {
    result.valid = false;
    // Get the detailed error message from SFCGAL
    result.reason = validity.reason();

    // Add additional context
    std::ostringstream details;
    details << "Geometry type: " << geom.geometryType() << "\n";
    details << "WKT: " << geom.asText(5) << "\n";

    if (geom.isEmpty()) {
      details << "Additional info: Geometry is empty\n";
    }

    if (geom.is3D()) {
      details << "Additional info: 3D geometry\n";
    }

    if (geom.isMeasured()) {
      details << "Additional info: Has M coordinates\n";
    }

    result.details = details.str();
  }

  return result;
}

auto
validate_operation_inputs(const std::string      &operation,
                          const SFCGAL::Geometry *first_geom,
                          const SFCGAL::Geometry *second_geom)
    -> ValidationResult
{
  ValidationResult result{true, "", ""};

  if (first_geom == nullptr) {
    result.valid  = false;
    result.reason = "First geometry is null";
    return result;
  }

  // Validate first geometry
  auto val_a = validate_geometry(*first_geom);
  if (!val_a.valid) {
    result.valid   = false;
    result.reason  = "First geometry invalid: " + val_a.reason;
    result.details = val_a.details;
    return result;
  }

  // Check if operation requires second geometry using registry lookup
  if (Operations::operation_requires_second_geometry(operation)) {
    if (second_geom == nullptr) {
      result.valid  = false;
      result.reason = "Operation '" + operation + "' requires second geometry";
      return result;
    }

    auto val_b = validate_geometry(*second_geom);
    if (!val_b.valid) {
      result.valid   = false;
      result.reason  = "Second geometry invalid: " + val_b.reason;
      result.details = val_b.details;
    }
  }

  return result;
}

void
report_error(const std::string &context, const std::exception &exception)
{
  std::cerr << "\033[1;31mError\033[0m";
  if (!context.empty()) {
    std::cerr << " in " << context;
  }
  std::cerr << ": " << exception.what() << "\n";
}

void
report_validation_error(const ValidationResult &result)
{
  std::cerr << "\033[1;31mValidation Error:\033[0m " << result.reason << "\n";
  if (!result.details.empty()) {
    std::cerr << "\033[33mDetails:\033[0m\n" << result.details;
  }
}

void
print_geometry_info(const SFCGAL::Geometry &geom)
{
  std::cout << "\033[1;36mGeometry Information:\033[0m\n";
  std::cout << "  Type: " << geom.geometryType() << "\n";
  std::cout << "  Dimension: " << geom.dimension() << "D\n";
  std::cout << "  Is 3D: " << (geom.is3D() ? "Yes" : "No") << "\n";
  std::cout << "  Is measured: " << (geom.isMeasured() ? "Yes" : "No") << "\n";
  std::cout << "  Is empty: " << (geom.isEmpty() ? "Yes" : "No") << "\n";
  std::cout << "  Is simple: "
            << (SFCGAL::algorithm::isSimple(geom) ? "Yes" : "No") << "\n";
  std::cout << "  Is valid: "
            << (SFCGAL::algorithm::isValid(geom) ? "Yes" : "No") << "\n";

  // Bounding box
  try {
    SFCGAL::Envelope envelope = geom.envelope();
    auto             poly     = envelope.toPolygon();
    if (poly) {
      std::cout << "  Bounding box: " << poly->asText(2) << "\n";
    }
  } catch (...) {
    // Envelope might not be available for all geometries
    std::cerr << "  Bounding box: Could not compute envelope\n";
  }
}

} // namespace ErrorHandler
