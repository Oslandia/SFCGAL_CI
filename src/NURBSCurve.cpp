// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/Envelope.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryVisitor.h"
#include "SFCGAL/algorithm/distance.h"
#include "SFCGAL/algorithm/offset.h"
#include "SFCGAL/detail/EnvelopeVisitor.h"
#include "SFCGAL/detail/ublas.h"
#include <CGAL/Bbox_3.h>
#include <algorithm>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <cmath>
#include <map>
#include <numeric>
#include <utility>

namespace SFCGAL {

/*
 * NURBS Curve Implementation Notes:
 *
 * This implementation follows algorithms from "The NURBS Book" by Piegl &
 * Tiller and provides full support for rational B-splines with XYZM
 * coordinates.
 *
 * Key Features:
 * - Exact evaluation using CGAL's exact arithmetic types
 * - Support for all coordinate types (XY, XYZ, XYM, XYZM)
 * - Rational weights handled in homogeneous space
 * - Multiple end conditions for curve fitting
 * - Parameter generation methods for stability
 *
 */

// Type aliases for convenience in static helper functions
using FT = NURBSCurve::FT; ///< Floating-point type alias for NURBS calculations

NURBSCurve::NURBSCurve() : _fitTolerance(FT(0)) {}

NURBSCurve::NURBSCurve(const std::vector<Point> &controlPoints,
                       unsigned int degree, KnotMethod knotMethod)
    : _controlPoints(controlPoints), _degree(degree), _fitTolerance(FT(0))
{
  if (!_controlPoints.empty()) {
    _weights.resize(_controlPoints.size(), FT(1));
    _knotVector = generateKnotVector(_controlPoints, _degree, knotMethod);

    auto [isValid, reason] = validateData();
    if (!isValid) {
      BOOST_THROW_EXCEPTION(Exception("Invalid NURBS curve data: " + reason));
    }
  }
}

NURBSCurve::NURBSCurve(const std::vector<Point> &controlPoints,
                       const std::vector<FT> &weights, unsigned int degree,
                       KnotMethod knotMethod)
    : _controlPoints(controlPoints), _weights(weights), _degree(degree),
      _fitTolerance(FT(0))
{
  if (!_controlPoints.empty()) {
    _knotVector = generateKnotVector(_controlPoints, _degree, knotMethod);

    auto [isValid, reason] = validateData();
    if (!isValid) {
      BOOST_THROW_EXCEPTION(Exception("Invalid NURBS curve data: " + reason));
    }
  }
}

NURBSCurve::NURBSCurve(const std::vector<Point> &controlPoints,
                       const std::vector<FT> &weights, unsigned int degree,
                       const std::vector<Knot> &knotVector)
    : _controlPoints(controlPoints), _weights(weights), _degree(degree),
      _knotVector(knotVector), _fitTolerance(FT(0))
{
  auto [isValid, reason] = validateData();
  if (!isValid) {
    BOOST_THROW_EXCEPTION(Exception("Invalid NURBS curve data: " + reason));
  }
}

auto
NURBSCurve::fromBezier(const std::vector<Point> &controlPoints)
    -> std::unique_ptr<NURBSCurve>
{
  if (controlPoints.empty()) {
    return std::make_unique<NURBSCurve>();
  }

  auto degree = static_cast<unsigned int>(controlPoints.size() - 1);

  std::vector<Knot> knots;
  knots.reserve(static_cast<size_t>(2) * (degree + 1));

  for (unsigned int idx = 0; idx <= degree; ++idx) {
    knots.emplace_back(0);
  }
  for (unsigned int idx = 0; idx <= degree; ++idx) {
    knots.emplace_back(1);
  }

  std::vector<FT> weights(controlPoints.size(), FT(1));
  return std::make_unique<NURBSCurve>(controlPoints, weights, degree, knots);
}

auto
NURBSCurve::fromBSpline(const std::vector<Point> &controlPoints,
                        unsigned int degree, const std::vector<Knot> &knots)
    -> std::unique_ptr<NURBSCurve>
{
  std::vector<FT> weights(controlPoints.size(), FT(1));
  return std::make_unique<NURBSCurve>(controlPoints, weights, degree, knots);
}

auto
NURBSCurve::createBSpline(const std::vector<Point> &controlPoints,
                          unsigned int degree) -> std::unique_ptr<NURBSCurve>
{
  if (controlPoints.empty()) {
    return std::make_unique<NURBSCurve>();
  }

  auto knots = generateKnotVector(controlPoints, degree, KnotMethod::UNIFORM);
  return fromBSpline(controlPoints, degree, knots);
}

// NOLINTBEGIN(readability-function-cognitive-complexity)
auto
NURBSCurve::createCircularArc(const Point &center, const FT &radius,
                              const FT &startAngle, const FT &endAngle,
                              const Point & /*normal*/,
                              KnotMethod /*knotMethod*/)
    -> std::unique_ptr<NURBSCurve>
{
  if (radius <= FT(0)) {
    BOOST_THROW_EXCEPTION(Exception("Arc radius must be positive"));
  }

  double startAngleDouble = CGAL::to_double(startAngle);
  double endAngleDouble   = CGAL::to_double(endAngle);
  double radiusDouble     = CGAL::to_double(radius);

  // Normalize angle range
  double angleSpan = endAngleDouble - startAngleDouble;
  while (angleSpan < 0.0) {
    angleSpan += 2.0 * M_PI;
  }
  while (angleSpan > 2.0 * M_PI) {
    angleSpan -= 2.0 * M_PI;
  }

  CoordinateType dimType = COORDINATE_XY;
  if (center.is3D() && center.isMeasured()) {
    dimType = COORDINATE_XYZM;
  } else if (center.is3D()) {
    dimType = COORDINATE_XYZ;
  } else if (center.isMeasured()) {
    dimType = COORDINATE_XYM;
  }

  std::vector<Point> controlPoints;
  std::vector<FT>    weights;

  // For arcs <= π, use standard quadratic rational representation
  // Handle half-circle degeneracy: when angleSpan ≈ π, cos(π/2) ≈ 0 creates
  // degenerate weight
  if (angleSpan <= M_PI + EPSILON) {
    // Check for near-π case (half-circle degeneracy)
    if (std::abs(angleSpan - M_PI) < EPSILON) {
      // For exact half-circle, use two quarter-circle segments to avoid
      // degeneracy
      double quarterAngle      = startAngleDouble + (angleSpan / 4.0);
      double midAngle          = startAngleDouble + (angleSpan / 2.0);
      double threeQuarterAngle = startAngleDouble + (3.0 * angleSpan / 4.0);

      // First quarter segment
      double cosQuarter = std::cos(M_PI / 8.0); // cos(π/8) for quarter circle

      // Control points for first quarter
      controlPoints.emplace_back(
          center.x() + FT(radiusDouble * std::cos(startAngleDouble)),
          center.y() + FT(radiusDouble * std::sin(startAngleDouble)),
          center.is3D() ? center.z() : FT(0),
          center.isMeasured() ? center.m() : NaN(), dimType);
      weights.emplace_back(1.0);

      controlPoints.emplace_back(
          center.x() + FT(radiusDouble * std::cos(quarterAngle) / cosQuarter),
          center.y() + FT(radiusDouble * std::sin(quarterAngle) / cosQuarter),
          center.is3D() ? center.z() : FT(0),
          center.isMeasured() ? center.m() : NaN(), dimType);
      weights.emplace_back(cosQuarter);

      controlPoints.emplace_back(
          center.x() + FT(radiusDouble * std::cos(midAngle)),
          center.y() + FT(radiusDouble * std::sin(midAngle)),
          center.is3D() ? center.z() : FT(0),
          center.isMeasured() ? center.m() : NaN(), dimType);
      weights.emplace_back(1.0);

      // Second quarter segment (continue the curve)
      controlPoints.emplace_back(
          center.x() +
              FT(radiusDouble * std::cos(threeQuarterAngle) / cosQuarter),
          center.y() +
              FT(radiusDouble * std::sin(threeQuarterAngle) / cosQuarter),
          center.is3D() ? center.z() : FT(0),
          center.isMeasured() ? center.m() : NaN(), dimType);
      weights.emplace_back(cosQuarter);

      controlPoints.emplace_back(
          center.x() + FT(radiusDouble * std::cos(endAngleDouble)),
          center.y() + FT(radiusDouble * std::sin(endAngleDouble)),
          center.is3D() ? center.z() : FT(0),
          center.isMeasured() ? center.m() : NaN(), dimType);
      weights.emplace_back(1.0);

      // Knot vector for degree 2 curve with 5 control points
      std::vector<Knot> knots = {FT(0),   FT(0), FT(0), FT(0.5),
                                 FT(0.5), FT(1), FT(1), FT(1)};
      return std::make_unique<NURBSCurve>(controlPoints, weights, 2, knots);
    }

    double halfSpan    = angleSpan / 2.0;
    double midAngle    = (startAngleDouble + endAngleDouble) / 2.0;
    double cosHalfSpan = std::cos(halfSpan);

    // Start point
    controlPoints.emplace_back(
        center.x() + FT(radiusDouble * std::cos(startAngleDouble)),
        center.y() + FT(radiusDouble * std::sin(startAngleDouble)),
        center.is3D() ? center.z() : FT(0),
        center.isMeasured() ? center.m() : NaN(), dimType);
    weights.emplace_back(1.0);

    // Middle control point (off the circle)
    controlPoints.emplace_back(
        center.x() + FT(radiusDouble * std::cos(midAngle) / cosHalfSpan),
        center.y() + FT(radiusDouble * std::sin(midAngle) / cosHalfSpan),
        center.is3D() ? center.z() : FT(0),
        center.isMeasured() ? center.m() : NaN(), dimType);
    weights.emplace_back(cosHalfSpan);

    // End point
    controlPoints.emplace_back(
        center.x() + FT(radiusDouble * std::cos(endAngleDouble)),
        center.y() + FT(radiusDouble * std::sin(endAngleDouble)),
        center.is3D() ? center.z() : FT(0),
        center.isMeasured() ? center.m() : NaN(), dimType);
    weights.emplace_back(1.0);

    std::vector<Knot> knots = {FT(0), FT(0), FT(0), FT(1), FT(1), FT(1)};
    return std::make_unique<NURBSCurve>(controlPoints, weights, 2, knots);
  }

  // For arcs > π, split into two segments
  double midAngle = startAngleDouble + (angleSpan / 2.0);

  // First segment: start to mid
  double halfSpan1    = (midAngle - startAngleDouble) / 2.0;
  double midAngle1    = (startAngleDouble + midAngle) / 2.0;
  double cosHalfSpan1 = std::cos(halfSpan1);

  controlPoints.emplace_back(
      center.x() + FT(radiusDouble * std::cos(startAngleDouble)),
      center.y() + FT(radiusDouble * std::sin(startAngleDouble)),
      center.is3D() ? center.z() : FT(0),
      center.isMeasured() ? center.m() : NaN(), dimType);
  weights.emplace_back(1.0);

  controlPoints.emplace_back(
      center.x() + FT(radiusDouble * std::cos(midAngle1) / cosHalfSpan1),
      center.y() + FT(radiusDouble * std::sin(midAngle1) / cosHalfSpan1),
      center.is3D() ? center.z() : FT(0),
      center.isMeasured() ? center.m() : NaN(), dimType);
  weights.emplace_back(cosHalfSpan1);

  controlPoints.emplace_back(center.x() + FT(radiusDouble * std::cos(midAngle)),
                             center.y() + FT(radiusDouble * std::sin(midAngle)),
                             center.is3D() ? center.z() : FT(0),
                             center.isMeasured() ? center.m() : NaN(), dimType);
  weights.emplace_back(1.0);

  // Second segment: mid to end
  double halfSpan2    = (endAngleDouble - midAngle) / 2.0;
  double midAngle2    = (midAngle + endAngleDouble) / 2.0;
  double cosHalfSpan2 = std::cos(halfSpan2);

  controlPoints.emplace_back(
      center.x() + FT(radiusDouble * std::cos(midAngle2) / cosHalfSpan2),
      center.y() + FT(radiusDouble * std::sin(midAngle2) / cosHalfSpan2),
      center.is3D() ? center.z() : FT(0),
      center.isMeasured() ? center.m() : NaN(), dimType);
  weights.emplace_back(cosHalfSpan2);

  controlPoints.emplace_back(
      center.x() + FT(radiusDouble * std::cos(endAngleDouble)),
      center.y() + FT(radiusDouble * std::sin(endAngleDouble)),
      center.is3D() ? center.z() : FT(0),
      center.isMeasured() ? center.m() : NaN(), dimType);
  weights.emplace_back(1.0);

  // Knot vector for two Bézier segments
  std::vector<Knot> knots = {FT(0),   FT(0), FT(0), FT(0.5),
                             FT(0.5), FT(1), FT(1), FT(1)};

  return std::make_unique<NURBSCurve>(controlPoints, weights, 2, knots);
}
// NOLINTEND(readability-function-cognitive-complexity)

// Helper functions for different end condition interpolations

// Static helper functions for matrix-based interpolation
static auto
findKnotSpanImpl(const FT &parameter, unsigned int degree,
                 const std::vector<FT> &knots, size_t numControlPoints)
    -> size_t
{
  // Handle boundary cases - parameter at or beyond curve end
  if (parameter >= knots[numControlPoints]) {
    return numControlPoints - 1;
  }

  // Handle boundary cases - parameter at or before curve start
  if (parameter <= knots[degree]) {
    return degree;
  }

  // Binary search for the knot span
  size_t low  = degree;
  size_t high = numControlPoints;
  size_t mid  = (low + high) / 2;

  while (parameter < knots[mid] || parameter >= knots[mid + 1]) {
    if (parameter < knots[mid]) {
      high = mid;
    } else {
      low = mid;
    }
    mid = (low + high) / 2;
  }

  return mid;
}

static auto
findKnotSpan(const FT &parameter, unsigned int degree,
             const std::vector<FT> &knots) -> size_t
{
  size_t numControlPoints = knots.size() - degree - 1;
  return findKnotSpanImpl(parameter, degree, knots, numControlPoints);
}

/**
 * Compute NURBS basis functions using Cox-de Boor recursion formula.
 *
 * From Piegl & Tiller Algorithm A2.2:
 * N_{i,0}(u) = 1 if u_i ≤ u < u_{i+1}, 0 otherwise
 * N_{i,p}(u) = [(u-u_i)/(u_{i+p}-u_i)] * N_{i,p-1}(u)
 *            + [(u_{i+p+1}-u)/(u_{i+p+1}-u_{i+1})] * N_{i+1,p-1}(u)
 *
 * Returns degree+1 non-zero basis function values at given parameter.
 */
static auto
computeBasisFunctions(size_t span, const FT &parameter, unsigned int degree,
                      const std::vector<FT> &knots) -> std::vector<FT>
{
  std::vector<FT> basis(degree + 1);
  std::vector<FT> left(degree + 1);
  std::vector<FT> right(degree + 1);

  basis[0] = FT(1.0);

  for (unsigned int j = 1; j <= degree; ++j) {
    left[j]  = parameter - knots[span + 1 - j];
    right[j] = knots[span + j] - parameter;

    FT saved = FT(0.0);
    for (unsigned int rIndex = 0; rIndex < j; ++rIndex) {
      FT temp       = basis[rIndex] / (right[rIndex + 1] + left[j - rIndex]);
      basis[rIndex] = saved + right[rIndex + 1] * temp;
      saved         = left[j - rIndex] * temp;
    }
    basis[j] = saved;
  }

  return basis;
}

// NOLINTBEGIN(readability-function-cognitive-complexity)
static auto
computeBasisDerivatives(size_t span, const FT &parameter, unsigned int degree,
                        const std::vector<FT> &knots,
                        unsigned int           maxDerivative)
    -> std::vector<std::vector<FT>>
{
  std::vector<std::vector<FT>> ders(maxDerivative + 1,
                                    std::vector<FT>(degree + 1, FT(0)));
  std::vector<FT>              left(degree + 1);
  std::vector<FT>              right(degree + 1);
  std::vector<std::vector<FT>> ndu(degree + 1, std::vector<FT>(degree + 1));
  std::vector<std::vector<FT>> a(2, std::vector<FT>(degree + 1));

  ndu[0][0] = FT(1.0);

  for (unsigned int j = 1; j <= degree; ++j) {
    left[j]  = parameter - knots[span + 1 - j];
    right[j] = knots[span + j] - parameter;
    FT saved = FT(0.0);

    for (unsigned int rIndex = 0; rIndex < j; ++rIndex) {
      ndu[j][rIndex] = right[rIndex + 1] + left[j - rIndex];
      FT temp = (ndu[j][rIndex] != FT(0)) ? ndu[rIndex][j - 1] / ndu[j][rIndex]
                                          : FT(0);
      ndu[rIndex][j] = saved + right[rIndex + 1] * temp;
      saved          = left[j - rIndex] * temp;
    }

    ndu[j][j] = saved;
  }

  for (unsigned int j = 0; j <= degree; ++j) {
    ders[0][j] = ndu[j][degree];
  }

  for (int rIdx = 0; rIdx <= static_cast<int>(degree); ++rIdx) {
    int s1Index = 0;
    int s2Index = 1;
    a[0][0]     = FT(1.0);

    for (unsigned int k = 1; k <= maxDerivative; ++k) {
      FT  derivative = FT(0.0);
      int rkValue    = rIdx - static_cast<int>(k);
      int pkValue    = static_cast<int>(degree) - static_cast<int>(k);

      if (rIdx >= static_cast<int>(k)) {
        a[s2Index][0] = (ndu[pkValue + 1][rkValue] != FT(0))
                            ? a[s1Index][0] / ndu[pkValue + 1][rkValue]
                            : FT(0);
        derivative    = a[s2Index][0] * ndu[rkValue][pkValue];
      }

      int j1Value = (rkValue >= -1) ? 1 : -rkValue;
      int j2Value = (rIdx - 1 <= pkValue) ? static_cast<int>(k) - 1
                                          : static_cast<int>(degree) - rIdx;

      for (int jIdx = j1Value; jIdx <= j2Value; ++jIdx) {
        a[s2Index][jIdx] = (ndu[pkValue + 1][rkValue + jIdx] != FT(0))
                               ? (a[s1Index][jIdx] - a[s1Index][jIdx - 1]) /
                                     ndu[pkValue + 1][rkValue + jIdx]
                               : FT(0);
        derivative += a[s2Index][jIdx] * ndu[rkValue + jIdx][pkValue];
      }

      if (rIdx <= pkValue) {
        a[s2Index][k] = (ndu[pkValue + 1][rIdx] != FT(0))
                            ? -a[s1Index][k - 1] / ndu[pkValue + 1][rIdx]
                            : FT(0);
        derivative += a[s2Index][k] * ndu[rIdx][pkValue];
      }

      ders[k][rIdx] = derivative;

      // Switch rows
      std::swap(s1Index, s2Index);
    }
  }

  // Multiply through by the correct factors
  int multiplier = static_cast<int>(degree);
  for (unsigned int k = 1; k <= maxDerivative; ++k) {
    for (unsigned int j = 0; j <= degree; ++j) {
      ders[k][j] *= FT(multiplier);
    }
    multiplier *= (static_cast<int>(degree) - static_cast<int>(k));
  }

  return ders;
}
// NOLINTEND(readability-function-cognitive-complexity)

auto
NURBSCurve::basisFunctionDerivatives(size_t span, const Parameter &parameter,
                                     unsigned int maxDerivative) const
    -> std::vector<std::vector<FT>>
{
  // Delegate to the existing free helper function
  return computeBasisDerivatives(span, parameter, _degree, _knotVector,
                                 maxDerivative);
}

auto
NURBSCurve::generateKnotVectorForEndCondition(
    const std::vector<Parameter> &parameters, unsigned int degree,
    EndCondition endCondition) -> std::vector<Knot>
{
  std::vector<Knot> knots;
  size_t            numKnots = parameters.size() + degree + 1;
  knots.reserve(numKnots);

  switch (endCondition) {
  case EndCondition::CLAMPED:
  case EndCondition::NATURAL:
  case EndCondition::TANGENT:
    // Standard clamped knot vector
    for (unsigned int idx = 0; idx <= degree; ++idx) {
      knots.push_back(parameters.front());
    }

    for (size_t pointIdx = 1; pointIdx < parameters.size() - degree;
         ++pointIdx) {
      FT sum = FT(0);
      for (unsigned int degIdx = 0; degIdx < degree; ++degIdx) {
        sum += parameters[pointIdx + degIdx];
      }
      knots.push_back(sum / FT(degree));
    }

    for (unsigned int idx = 0; idx <= degree; ++idx) {
      knots.push_back(parameters.back());
    }
    break;

  case EndCondition::PERIODIC:
    // Handled separately in interpolatePeriodicCurve
    BOOST_THROW_EXCEPTION(Exception("Periodic knots handled separately"));
    break;
  }

  return knots;
}

auto
NURBSCurve::interpolateClampedCurve(const std::vector<Point>     &points,
                                    const std::vector<Parameter> &parameters,
                                    unsigned int                  degree,
                                    const std::vector<Knot>      &knots)
    -> std::vector<Point>
{
  using namespace detail::ublas;

  size_t numPoints = points.size();
  if (numPoints < degree + 1) {
    BOOST_THROW_EXCEPTION(
        Exception("Not enough points for clamped interpolation"));
  }

  // For clamped interpolation, we solve the system N * P = Q
  // where N is the basis function matrix, P are control points, Q are data
  // points The clamped condition ensures the curve passes exactly through
  // endpoints

  matrix<double> basisMatrix(numPoints, numPoints);

  // Fill basis function matrix
  for (size_t i = 0; i < numPoints; ++i) {
    const FT &parameter = parameters[i];

    // Find knot span
    size_t span = findKnotSpan(parameter, degree, knots);

    // Compute basis functions at this parameter
    auto basis = computeBasisFunctions(span, parameter, degree, knots);

    // Initialize row to zero
    std::fill(basisMatrix.data().begin() + (i * numPoints),
              basisMatrix.data().begin() + ((i + 1) * numPoints), 0.0);

    // Fill non-zero entries in this row
    size_t baseIdx = span - degree;
    for (unsigned int j = 0; j <= degree && baseIdx + j < numPoints; ++j) {
      basisMatrix(i, baseIdx + j) = CGAL::to_double(basis[j]);
    }
  }

  // Solve the system for each coordinate dimension
  std::vector<Point> controlPoints;
  controlPoints.reserve(numPoints);

  // Determine coordinate dimensions
  bool is3D       = points.front().is3D();
  bool isMeasured = points.front().isMeasured();

  try {
    // Solve system using LU decomposition for each coordinate
    permutation_matrix<std::size_t> pm(numPoints);
    matrix<double>                  NCopy;

    // Solve for X coordinates
    vector<double> qx(numPoints);
    for (size_t i = 0; i < numPoints; ++i) {
      qx(i) = CGAL::to_double(points[i].x());
    }

    NCopy = basisMatrix;
    lu_factorize(NCopy, pm);
    lu_substitute(NCopy, pm, qx);
    vector<double> px = qx;

    // Solve for Y coordinates
    vector<double> qy(numPoints);
    for (size_t i = 0; i < numPoints; ++i) {
      qy(i) = CGAL::to_double(points[i].y());
    }

    NCopy = basisMatrix;
    lu_factorize(NCopy, pm);
    lu_substitute(NCopy, pm, qy);
    vector<double> py = qy;

    // Solve for Z coordinates if 3D
    vector<double> pz(numPoints);
    if (is3D) {
      vector<double> qz(numPoints);
      for (size_t i = 0; i < numPoints; ++i) {
        qz(i) = CGAL::to_double(points[i].z());
      }

      NCopy = basisMatrix;
      lu_factorize(NCopy, pm);
      lu_substitute(NCopy, pm, qz);
      pz = qz;
    }

    // Solve for M coordinates if measured
    vector<double> pm_coord(numPoints);
    if (isMeasured) {
      vector<double> qm(numPoints);
      for (size_t i = 0; i < numPoints; ++i) {
        qm(i) = CGAL::to_double(points[i].m());
      }

      NCopy = basisMatrix;
      lu_factorize(NCopy, pm);
      lu_substitute(NCopy, pm, qm);
      pm_coord = qm;
    }

    // Build control points
    for (size_t i = 0; i < numPoints; ++i) {
      FT x = FT(px(i));
      FT y = FT(py(i));

      if (is3D && isMeasured) {
        FT     z = FT(pz(i));
        double m = CGAL::to_double(FT(pm_coord(i)));
        controlPoints.emplace_back(x, y, z,
                                   m); // Uses (FT, FT, FT, double) constructor
      } else if (is3D) {
        FT z = FT(pz(i));
        controlPoints.emplace_back(x, y, z);
      } else if (isMeasured) {
        double x_d = CGAL::to_double(x);
        double y_d = CGAL::to_double(y);
        double m   = pm_coord(i);
        controlPoints.emplace_back(
            x_d, y_d, 0.0, m,
            COORDINATE_XYM); // Use explicit constructor with CoordinateType
      } else {
        controlPoints.emplace_back(x, y);
      }
    }

  } catch (const std::exception &e) {
    BOOST_THROW_EXCEPTION(
        Exception("Failed to solve clamped interpolation system: " +
                  std::string(e.what())));
  }

  return controlPoints;
}

// NOLINTBEGIN(readability-function-cognitive-complexity)
auto
NURBSCurve::interpolateNaturalCurve(const std::vector<Point>     &points,
                                    const std::vector<Parameter> &parameters,
                                    unsigned int                  degree,
                                    const std::vector<Knot>      &knots)
    -> std::vector<Point>
{
  using namespace detail::ublas;

  size_t numDataPoints = points.size();
  if (numDataPoints < degree + 1) {
    BOOST_THROW_EXCEPTION(
        Exception("Not enough points for natural interpolation"));
  }

  // Set up the interpolation matrix system N * P = Q
  // where N is the basis function matrix, P are control points, Q are data
  // points

  matrix<double> naturalBasisMatrix(numDataPoints, numDataPoints);

  // Fill basis function matrix
  for (size_t i = 0; i < numDataPoints; ++i) {
    const FT &parameter = parameters[i];

    // Find knot span
    size_t span = findKnotSpan(parameter, degree, knots);

    // Compute basis functions
    auto basis = computeBasisFunctions(span, parameter, degree, knots);

    // Fill row of matrix
    std::fill(naturalBasisMatrix.data().begin() + (i * numDataPoints),
              naturalBasisMatrix.data().begin() + ((i + 1) * numDataPoints),
              0.0);

    size_t baseIdx = span - degree;
    for (unsigned int j = 0; j <= degree && baseIdx + j < numDataPoints; ++j) {
      naturalBasisMatrix(i, baseIdx + j) = CGAL::to_double(basis[j]);
    }
  }

  // Add natural end conditions (minimize curvature at ends)
  if (degree >= 2) {
    // Replace first and last equations with curvature minimization
    // Second derivative at ends should be minimal

    // For simplicity, use zero second derivative conditions
    // This is a simplified natural condition - full implementation would
    // minimize integral of curvature

    std::fill(naturalBasisMatrix.data().begin(),
              naturalBasisMatrix.data().begin() + numDataPoints, 0.0);
    std::fill(naturalBasisMatrix.data().begin() +
                  ((numDataPoints - 1) * numDataPoints),
              naturalBasisMatrix.data().begin() +
                  (numDataPoints * numDataPoints),
              0.0);

    // Set up second derivative = 0 at start
    const FT &startParam = parameters.front();
    size_t    startSpan  = findKnotSpan(startParam, degree, knots);
    auto      startBasis2nd =
        computeBasisDerivatives(startSpan, startParam, degree, knots, 2);

    size_t startBaseIdx = startSpan - degree;
    for (unsigned int j = 0; j <= degree && startBaseIdx + j < numDataPoints;
         ++j) {
      naturalBasisMatrix(0, startBaseIdx + j) =
          CGAL::to_double(startBasis2nd[2][j]);
    }

    // Set up second derivative = 0 at end
    const FT &endParam = parameters.back();
    size_t    endSpan  = findKnotSpan(endParam, degree, knots);
    auto      endBasis2nd =
        computeBasisDerivatives(endSpan, endParam, degree, knots, 2);

    size_t endBaseIdx = endSpan - degree;
    for (unsigned int j = 0; j <= degree && endBaseIdx + j < numDataPoints;
         ++j) {
      naturalBasisMatrix(numDataPoints - 1, endBaseIdx + j) =
          CGAL::to_double(endBasis2nd[2][j]);
    }
  }

  // Solve for each coordinate dimension
  std::vector<Point> controlPoints;
  controlPoints.reserve(numDataPoints);

  // Determine coordinate type from first point
  CoordinateType coordType = COORDINATE_XY;
  if (points[0].is3D() && points[0].isMeasured()) {
    coordType = COORDINATE_XYZM;
  } else if (points[0].is3D()) {
    coordType = COORDINATE_XYZ;
  } else if (points[0].isMeasured()) {
    coordType = COORDINATE_XYM;
  }

  // Solve for X coordinates
  vector<double> qx(numDataPoints);
  vector<double> px(numDataPoints);
  for (size_t i = 0; i < numDataPoints; ++i) {
    qx(i) = CGAL::to_double(points[i].x());
    if (degree >= 2 && (i == 0 || i == numDataPoints - 1)) {
      qx(i) = 0.0; // Natural end condition: zero second derivative
    }
  }

  // Solve N * px = qx
  permutation_matrix<size_t> pm(numDataPoints);
  matrix<double>             NCopy = naturalBasisMatrix;
  lu_factorize(NCopy, pm);
  lu_substitute(NCopy, pm, qx);
  px = qx; // qx now contains solution

  // Solve for Y coordinates
  vector<double> qy(numDataPoints);
  vector<double> py(numDataPoints);
  for (size_t i = 0; i < numDataPoints; ++i) {
    qy(i) = CGAL::to_double(points[i].y());
    if (degree >= 2 && (i == 0 || i == numDataPoints - 1)) {
      qy(i) = 0.0; // Natural end condition
    }
  }

  NCopy = naturalBasisMatrix;
  lu_factorize(NCopy, pm);
  lu_substitute(NCopy, pm, qy);
  py = qy;

  // Solve for Z coordinates if 3D
  vector<double> pz(numDataPoints);
  if (points[0].is3D()) {
    vector<double> qz(numDataPoints);
    for (size_t i = 0; i < numDataPoints; ++i) {
      qz(i) = CGAL::to_double(points[i].z());
      if (degree >= 2 && (i == 0 || i == numDataPoints - 1)) {
        qz(i) = 0.0;
      }
    }

    NCopy = naturalBasisMatrix;
    lu_factorize(NCopy, pm);
    lu_substitute(NCopy, pm, qz);
    pz = qz;
  }

  // Solve for M coordinates if measured
  vector<double> pm_coord(numDataPoints);
  if (points[0].isMeasured()) {
    vector<double> qm(numDataPoints);
    for (size_t i = 0; i < numDataPoints; ++i) {
      qm(i) = CGAL::to_double(points[i].m());
      if (degree >= 2 && (i == 0 || i == numDataPoints - 1)) {
        qm(i) = 0.0;
      }
    }

    NCopy = naturalBasisMatrix;
    lu_factorize(NCopy, pm);
    lu_substitute(NCopy, pm, qm);
    pm_coord = qm;
  }

  // Build control points
  for (size_t i = 0; i < numDataPoints; ++i) {
    FT     x = FT(px(i));
    FT     y = FT(py(i));
    FT     z = points[0].is3D() ? FT(pz(i)) : FT(0);
    double m = points[0].isMeasured() ? pm_coord(i) : NaN();

    controlPoints.emplace_back(x, y, z, m, coordType);
  }

  return controlPoints;
}
// NOLINTEND(readability-function-cognitive-complexity)

auto
NURBSCurve::interpolatePeriodicCurve(const std::vector<Point> &points,
                                     unsigned int              degree,
                                     KnotMethod /*knotMethod*/)
    -> std::unique_ptr<NURBSCurve>
{
  if (points.size() < degree + 1) {
    BOOST_THROW_EXCEPTION(
        Exception("Not enough points for periodic interpolation"));
  }

  // Check if curve is actually closed
  if (algorithm::distance(points.front(), points.back()) > EPSILON) {
    BOOST_THROW_EXCEPTION(
        Exception("Points must form closed curve for periodic interpolation"));
  }

  // For periodic curves, we need to wrap points and create uniform knot spacing
  std::vector<Point> periodicPoints = points;

  // Remove duplicate end point for internal processing
  if (periodicPoints.size() > 1) {
    periodicPoints.pop_back();
  }

  size_t numPeriodicPoints = periodicPoints.size();

  // Create periodic parameter vector
  std::vector<Parameter> parameters;
  parameters.reserve(numPeriodicPoints);

  for (size_t i = 0; i < numPeriodicPoints; ++i) {
    parameters.push_back(FT(i) / FT(numPeriodicPoints));
  }

  // Create periodic knot vector
  std::vector<Knot> knots;
  size_t            numKnots = numPeriodicPoints + degree + 1;
  knots.reserve(numKnots);

  // Periodic knots wrap around
  for (size_t i = 0; i < numKnots; ++i) {
    int knotIndex = static_cast<int>(i) - static_cast<int>(degree);
    knots.push_back(FT(knotIndex) / FT(numPeriodicPoints));
  }

  // For now, use simple periodic control points (same as data points)
  // A full implementation would solve a periodic interpolation system
  std::vector<Point> controlPoints = periodicPoints;

  auto curve = std::make_unique<NURBSCurve>(
      controlPoints, std::vector<FT>(controlPoints.size(), FT(1)), degree,
      knots);

  curve->_fitPoints    = points;
  curve->_fitTolerance = FT(0);

  return curve;
}

auto
NURBSCurve::interpolateCurve(const std::vector<Point> &points,
                             unsigned int degree, KnotMethod knotMethod,
                             EndCondition endCondition)
    -> std::unique_ptr<NURBSCurve>
{
  if (points.size() < 2) {
    BOOST_THROW_EXCEPTION(
        Exception("Need at least 2 points for interpolation"));
  }

  // Handle periodic curves specially
  if (endCondition == EndCondition::PERIODIC) {
    return interpolatePeriodicCurve(points, degree, knotMethod);
  }

  if (degree >= points.size()) {
    degree = static_cast<unsigned int>(points.size() - 1);
  }

  // Compute parameter values
  auto parameters = computeParameters(points, knotMethod);

  // Generate appropriate knot vector based on end condition
  std::vector<Knot> knots =
      generateKnotVectorForEndCondition(parameters, degree, endCondition);

  std::vector<Point> controlPoints;

  // Handle different end conditions
  switch (endCondition) {
  case EndCondition::CLAMPED:
    controlPoints = interpolateClampedCurve(points, parameters, degree, knots);
    break;

  case EndCondition::NATURAL:
    controlPoints = interpolateNaturalCurve(points, parameters, degree, knots);
    break;

  case EndCondition::TANGENT:
    // For now, fall back to clamped if no tangent vectors specified
    // TODO: Add support for user-specified tangent vectors
    controlPoints = interpolateClampedCurve(points, parameters, degree, knots);
    break;

  case EndCondition::PERIODIC:
    // Already handled above
    BOOST_THROW_EXCEPTION(Exception("Periodic case should be handled earlier"));
    break;

  default:
    BOOST_THROW_EXCEPTION(Exception("Unknown EndCondition"));
  }

  auto curve = std::make_unique<NURBSCurve>(
      controlPoints, std::vector<FT>(controlPoints.size(), FT(1)), degree,
      knots);

  curve->_fitPoints    = points;
  curve->_fitTolerance = FT(0);

  return curve;
}

auto
NURBSCurve::generateApproximationKnotVector(
    const std::vector<Parameter> &parameters, unsigned int degree,
    size_t numControlPoints) -> std::vector<Knot>
{
  /**
   * Generate knot vector for NURBS approximation using geomdl method.
   *
   * This is compute_knot_vector2 from geomdl which uses Piegl & Tiller
   * Equation 9.68 from "The NURBS Book".
   *
   * This method ensures that every knot span contains at least one parameter
   * value, which produces better approximations than simple averaging.
   */
  size_t num_dpts = parameters.size(); // Number of data points
  size_t num_cpts = numControlPoints;  // Number of control points

  std::vector<Knot> knots;
  knots.reserve(num_cpts + degree + 1);

  // First degree+1 knots are 0 (start knot with multiplicity degree+1)
  for (unsigned int i = 0; i <= degree; ++i) {
    knots.emplace_back(0);
  }

  // Compute internal knots using geomdl method (Equation 9.68)
  // d = n / (m - p) where n = num_data_points, m = num_control_points - 1, p =
  // degree
  double deltaParam =
      static_cast<double>(num_dpts) / static_cast<double>(num_cpts - degree);

  // Find internal knots
  for (size_t j = 1; j < num_cpts - degree; ++j) {
    double jDelta     = static_cast<double>(j) * deltaParam;
    int    paramIndex = static_cast<int>(jDelta);
    double alpha      = jDelta - static_cast<double>(paramIndex);

    // Ensure we don't go out of bounds
    if (paramIndex >= static_cast<int>(parameters.size())) {
      paramIndex = static_cast<int>(parameters.size() - 1);
      alpha      = 0.0;
    }
    paramIndex = std::max(paramIndex, 1);

    // Compute knot as linear interpolation between parameters[paramIndex-1] and
    // parameters[paramIndex]
    FT temp_kv = (FT(1.0 - alpha) * parameters[paramIndex - 1]) +
                 (FT(alpha) * parameters[paramIndex]);
    knots.push_back(temp_kv);
  }

  // Last degree+1 knots are 1 (end knot with multiplicity degree+1)
  for (unsigned int i = 0; i <= degree; ++i) {
    knots.emplace_back(1);
  }

  return knots;
}

// NOLINTBEGIN(readability-function-cognitive-complexity)
auto
NURBSCurve::approximateCurve(const std::vector<Point> &points,
                             unsigned int degree, FT tolerance,
                             size_t maxControlPoints)
    -> std::unique_ptr<NURBSCurve>
{
  if (points.size() < 2) {
    BOOST_THROW_EXCEPTION(
        Exception("Need at least 2 points for approximation"));
  }

  if (degree < 1) {
    BOOST_THROW_EXCEPTION(Exception("Degree must be at least 1"));
  }

  // Ensure we have enough control points for the degree
  size_t minControlPoints = degree + 1;
  size_t numControlPoints = std::min(maxControlPoints, points.size());
  numControlPoints        = std::max(numControlPoints, minControlPoints);

  if (numControlPoints >= points.size()) {
    // If we need as many control points as data points, use interpolation
    // instead
    return interpolateCurve(points, degree, KnotMethod::CHORD_LENGTH,
                            EndCondition::CLAMPED);
  }

  /**
   * NURBS Curve Approximation - Piegl & Tiller Algorithm A9.1
   *
   * This implements true least-squares approximation, generating control points
   * that minimize the sum of squared distances to input data points.
   *
   * Mathematical Foundation:
   * - Given: n+1 data points Q_k, k = 0,...,n
   * - Find: m+1 control points P_j that minimize ||C(u_k) - Q_k||²
   * - Where: C(u) = Σ_{j=0}^m N_{j,p}(u) * P_j (NURBS curve)
   *
   * The least-squares problem becomes: N^T * N * P = N^T * Q
   * where N_{i,j} = N_{j,p}(u_i) are the basis function values
   */

  using namespace detail::ublas;

  size_t numDataMinusOne =
      points.size() - 1; // Number of data points - 1 (index range)
  size_t numControlMinusOne =
      numControlPoints - 1; // Number of control points - 1

  // Step 1: Compute parameter values for data points (chord length
  // parameterization)
  std::vector<Parameter> parameters =
      computeParameters(points, KnotMethod::CHORD_LENGTH);

  // Step 2: Generate knot vector for approximation
  std::vector<Knot> knots = generateApproximationKnotVector(
      parameters, degree, numControlMinusOne + 1);

  // Step 3: Set up the least-squares system based on approximation mode
  // Mode determines whether to fix endpoints (SMOOTH) or include all points
  // (FAITHFUL)

  std::vector<Point> controlPoints;
  controlPoints.reserve(numControlMinusOne + 1);

  // Determine coordinate dimensions
  bool is3D       = points.front().is3D();
  bool isMeasured = points.front().isMeasured();

  try {
    // Use geomdl-like behavior:
    // Fix endpoints but still minimize error to ALL data points

    // Initialize control points with fixed endpoints
    controlPoints.resize(numControlMinusOne + 1);
    controlPoints[0] = points[0]; // Fix first control point to first data point
    controlPoints[numControlMinusOne] =
        points[numDataMinusOne]; // Fix last control point to last data point

    if (numControlMinusOne > 1) {
      // Build basis function matrix for INTERIOR data points only (like geomdl)
      // This excludes the first and last data points which are already
      // satisfied by the fixed control points
      matrix<double> N_interior(numDataMinusOne - 1,
                                numControlMinusOne -
                                    1); // Interior data points (1 to n-1) x
                                        // interior control points (1 to m-1)

      // Fill basis function matrix for INTERIOR data points only
      for (size_t i = 1; i < numDataMinusOne;
           ++i) { // Loop from 1 to n-1 (excluding endpoints)
        const Parameter &t     = parameters[i];
        size_t           span  = findKnotSpan(t, degree, knots);
        auto             basis = computeBasisFunctions(span, t, degree, knots);

        // Initialize row to zero
        for (size_t j = 0; j < numControlMinusOne - 1; ++j) {
          N_interior(i - 1, j) = 0.0; // i-1 because matrix starts at 0
        }

        // Fill non-zero basis functions for interior control points only
        size_t baseIdx = span - degree;
        for (unsigned int k = 0;
             k <= degree && baseIdx + k <= numControlMinusOne; ++k) {
          if (baseIdx + k >= 1 &&
              baseIdx + k < numControlMinusOne) { // Only interior control
                                                  // points (1 to m-1)
            N_interior(i - 1, baseIdx + k - 1) = CGAL::to_double(basis[k]);
          }
        }
      }

      // Compute modified RHS for interior points: Rk = Qk - N0*P0 - Nm*Pm
      // where P0 and Pm are fixed to data endpoints
      vector<double> rk_x(numDataMinusOne - 1);
      vector<double> rk_y(numDataMinusOne - 1);
      vector<double> rk_z(numDataMinusOne - 1, 0.0);
      vector<double> rk_m(numDataMinusOne - 1, 0.0);

      for (size_t i = 1; i < numDataMinusOne;
           ++i) { // Only interior data points
        const Parameter &t     = parameters[i];
        size_t           span  = findKnotSpan(t, degree, knots);
        auto             basis = computeBasisFunctions(span, t, degree, knots);

        // Start with the interior data point
        rk_x(i - 1) = CGAL::to_double(points[i].x());
        rk_y(i - 1) = CGAL::to_double(points[i].y());
        if (is3D) {
          rk_z(i - 1) = CGAL::to_double(points[i].z());
        }
        if (isMeasured) {
          rk_m(i - 1) = CGAL::to_double(points[i].m());
        }

        // Subtract contribution from fixed endpoints (like geomdl)
        size_t baseIdx = span - degree;
        for (unsigned int k = 0;
             k <= degree && baseIdx + k <= numControlMinusOne; ++k) {
          if (baseIdx + k == 0) {
            // Contribution from first control point (fixed to points[0])
            double firstBasisValue = CGAL::to_double(basis[k]);
            rk_x(i - 1) -= firstBasisValue * CGAL::to_double(points[0].x());
            rk_y(i - 1) -= firstBasisValue * CGAL::to_double(points[0].y());
            if (is3D) {
              rk_z(i - 1) -= firstBasisValue * CGAL::to_double(points[0].z());
            }
            if (isMeasured) {
              rk_m(i - 1) -= firstBasisValue * CGAL::to_double(points[0].m());
            }
          } else if (baseIdx + k == numControlMinusOne) {
            // Contribution from last control point (fixed to points[n])
            double lastBasisValue = CGAL::to_double(basis[k]);
            rk_x(i - 1) -=
                lastBasisValue * CGAL::to_double(points[numDataMinusOne].x());
            rk_y(i - 1) -=
                lastBasisValue * CGAL::to_double(points[numDataMinusOne].y());
            if (is3D) {
              rk_z(i - 1) -=
                  lastBasisValue * CGAL::to_double(points[numDataMinusOne].z());
            }
            if (isMeasured) {
              rk_m(i - 1) -=
                  lastBasisValue * CGAL::to_double(points[numDataMinusOne].m());
            }
          }
        }
      }

      // Solve the reduced least-squares system: N_interior^T * N_interior * P =
      // N_interior^T * Rk
      matrix<double>                  NTN = prod(trans(N_interior), N_interior);
      permutation_matrix<std::size_t> pm(numControlMinusOne - 1);

      // Solve for X coordinates
      vector<double> ntqx     = prod(trans(N_interior), rk_x);
      matrix<double> NTN_copy = NTN;
      lu_factorize(NTN_copy, pm);
      lu_substitute(NTN_copy, pm, ntqx);

      // Solve for Y coordinates
      vector<double> ntqy = prod(trans(N_interior), rk_y);
      NTN_copy            = NTN;
      lu_factorize(NTN_copy, pm);
      lu_substitute(NTN_copy, pm, ntqy);

      // Solve for Z coordinates if 3D
      vector<double> pz_solved(numControlMinusOne - 1, 0.0);
      if (is3D) {
        vector<double> ntqz = prod(trans(N_interior), rk_z);
        NTN_copy            = NTN;
        lu_factorize(NTN_copy, pm);
        lu_substitute(NTN_copy, pm, ntqz);
        pz_solved = ntqz;
      }

      // Solve for M coordinates if measured
      vector<double> pm_solved(numControlMinusOne - 1, 0.0);
      if (isMeasured) {
        vector<double> ntqm = prod(trans(N_interior), rk_m);
        NTN_copy            = NTN;
        lu_factorize(NTN_copy, pm);
        lu_substitute(NTN_copy, pm, ntqm);
        pm_solved = ntqm;
      }

      // Fill interior control points
      for (size_t i = 1; i < numControlMinusOne; ++i) {
        FT x = FT(ntqx(i - 1));
        FT y = FT(ntqy(i - 1));

        if (is3D && isMeasured) {
          FT     z         = FT(pz_solved(i - 1));
          double m_val     = pm_solved(i - 1);
          controlPoints[i] = Point(x, y, z, m_val);
        } else if (is3D) {
          FT z             = FT(pz_solved(i - 1));
          controlPoints[i] = Point(x, y, z);
        } else if (isMeasured) {
          double x_d       = CGAL::to_double(x);
          double y_d       = CGAL::to_double(y);
          double m_val     = pm_solved(i - 1);
          controlPoints[i] = Point(x_d, y_d, 0.0, m_val, COORDINATE_XYM);
        } else {
          controlPoints[i] = Point(x, y);
        }
      }
    }

  } catch (const std::exception &e) {
    BOOST_THROW_EXCEPTION(
        Exception("Failed to solve least-squares approximation system: " +
                  std::string(e.what())));
  }

  // Create the approximating curve
  auto curve = std::make_unique<NURBSCurve>(
      controlPoints, std::vector<FT>(controlPoints.size(), FT(1)), degree,
      knots);

  // Store fit information
  curve->_fitPoints    = points;
  curve->_fitTolerance = std::move(tolerance);

  return curve;
}
// NOLINTEND(readability-function-cognitive-complexity)

auto
NURBSCurve::fitCurve(const std::vector<Point> &points, unsigned int degree,
                     FitMethod fitMethod, KnotMethod knotMethod,
                     EndCondition endCondition, FT tolerance,
                     size_t maxControlPoints) -> std::unique_ptr<NURBSCurve>
{
  if (points.empty()) {
    BOOST_THROW_EXCEPTION(Exception("Cannot fit curve with no points"));
  }

  switch (fitMethod) {
  case FitMethod::INTERPOLATE:
    // Use interpolation - pass through all points exactly
    return interpolateCurve(points, degree, knotMethod, endCondition);

  case FitMethod::APPROXIMATE:
    // Use approximation - smooth curve within tolerance
    return approximateCurve(points, degree, std::move(tolerance),
                            maxControlPoints);

  default:
    BOOST_THROW_EXCEPTION(Exception("Unknown FitMethod"));
  }
}

void
NURBSCurve::accept(GeometryVisitor &visitor)
{
  visitor.visit(*this);
}

void
NURBSCurve::accept(ConstGeometryVisitor &visitor) const
{
  visitor.visit(*this);
}

auto
NURBSCurve::geometryTypeId() const -> GeometryType
{
  return TYPE_NURBSCURVE;
}

auto
NURBSCurve::geometryType() const -> std::string
{
  return "NURBSCurve";
}

auto
NURBSCurve::clone() const -> NURBSCurve *
{
  return new NURBSCurve(*this);
}

auto
NURBSCurve::isEmpty() const -> bool
{
  return _controlPoints.empty();
}

auto
NURBSCurve::is3D() const -> bool
{
  return !_controlPoints.empty() && _controlPoints[0].is3D();
}

auto
NURBSCurve::isMeasured() const -> bool
{
  return !_controlPoints.empty() && _controlPoints[0].isMeasured();
}

auto
NURBSCurve::dropZ() -> bool
{
  bool hadZ = is3D();
  for (auto &point : _controlPoints) {
    point.dropZ();
  }
  return hadZ;
}

auto
NURBSCurve::dropM() -> bool
{
  bool hadM = isMeasured();
  for (auto &point : _controlPoints) {
    point.dropM();
  }
  return hadM;
}

void
NURBSCurve::swapXY()
{
  for (auto &point : _controlPoints) {
    point.swapXY();
  }
}

auto
NURBSCurve::evaluate(Parameter parameter) const -> Point
{
  if (_controlPoints.empty()) {
    BOOST_THROW_EXCEPTION(Exception("Cannot evaluate empty NURBS curve"));
  }

  auto bounds = parameterBounds();
  if (parameter < bounds.first || parameter > bounds.second) {
    BOOST_THROW_EXCEPTION(Exception("Parameter outside valid range"));
  }

  if (_controlPoints.size() == 1) {
    return _controlPoints[0];
  }

  size_t span = findSpan(parameter);
  return deBoorRational(span, parameter);
}

auto
NURBSCurve::derivative(Parameter parameter, unsigned int order) const -> Point
{
  if (order == 0) {
    return evaluate(parameter);
  }

  auto derivatives = computeDerivatives(parameter, order);
  if (order < derivatives.size()) {
    return derivatives[order];
  }

  if (is3D() && isMeasured()) {
    return {FT(0), FT(0), FT(0), 0.0};
  }
  if (is3D()) {
    return {FT(0), FT(0), FT(0)};
  }
  if (isMeasured()) {
    Point point{FT(0), FT(0)};
    point.setM(0.0);
    return point;
  }
  return {FT(0), FT(0)};
}

auto
NURBSCurve::tangent(Parameter parameter) const -> Point
{
  Point firstDerivative = derivative(parameter, 1);

  FT magnitude = std::sqrt(CGAL::to_double(
      firstDerivative.x() * firstDerivative.x() +
      firstDerivative.y() * firstDerivative.y() +
      (firstDerivative.is3D() ? firstDerivative.z() * firstDerivative.z()
                              : FT(0))));

  if (magnitude < EPSILON) {
    BOOST_THROW_EXCEPTION(
        Exception("Cannot compute tangent at singular point"));
  }

  return {firstDerivative.x() / magnitude, firstDerivative.y() / magnitude,
          firstDerivative.is3D() ? firstDerivative.z() / magnitude : FT(0)};
}

auto
NURBSCurve::normal(Parameter parameter) const -> Point
{
  if (is3D()) {
    BOOST_THROW_EXCEPTION(
        Exception("Normal vector not uniquely defined for 3D curves"));
  }

  Point tangentVec = tangent(parameter);
  return {-tangentVec.y(), tangentVec.x()};
}

auto
NURBSCurve::binormal(Parameter parameter) const -> Point
{
  if (!is3D()) {
    BOOST_THROW_EXCEPTION(Exception("Binormal only defined for 3D curves"));
  }

  Point firstDeriv  = derivative(parameter, 1);
  Point secondDeriv = derivative(parameter, 2);

  FT crossX =
      firstDeriv.y() * secondDeriv.z() - firstDeriv.z() * secondDeriv.y();
  FT crossY =
      firstDeriv.z() * secondDeriv.x() - firstDeriv.x() * secondDeriv.z();
  FT crossZ =
      firstDeriv.x() * secondDeriv.y() - firstDeriv.y() * secondDeriv.x();

  FT magnitude = std::sqrt(
      CGAL::to_double(crossX * crossX + crossY * crossY + crossZ * crossZ));

  if (magnitude < EPSILON) {
    BOOST_THROW_EXCEPTION(
        Exception("Cannot compute binormal at inflection point"));
  }

  return {crossX / magnitude, crossY / magnitude, crossZ / magnitude};
}

auto
NURBSCurve::curvature(Parameter parameter) const -> FT
{
  Point firstDeriv  = derivative(parameter, 1);
  Point secondDeriv = derivative(parameter, 2);

  if (is3D()) {
    FT crossX =
        firstDeriv.y() * secondDeriv.z() - firstDeriv.z() * secondDeriv.y();
    FT crossY =
        firstDeriv.z() * secondDeriv.x() - firstDeriv.x() * secondDeriv.z();
    FT crossZ =
        firstDeriv.x() * secondDeriv.y() - firstDeriv.y() * secondDeriv.x();

    FT crossMagnitude = std::sqrt(
        CGAL::to_double(crossX * crossX + crossY * crossY + crossZ * crossZ));
    FT firstMagnitude = std::sqrt(CGAL::to_double(
        firstDeriv.x() * firstDeriv.x() + firstDeriv.y() * firstDeriv.y() +
        firstDeriv.z() * firstDeriv.z()));

    if (firstMagnitude < EPSILON) {
      return {0};
    }

    return crossMagnitude / (firstMagnitude * firstMagnitude * firstMagnitude);
  }
  FT numerator =
      firstDeriv.x() * secondDeriv.y() - firstDeriv.y() * secondDeriv.x();
  FT denominator =
      firstDeriv.x() * firstDeriv.x() + firstDeriv.y() * firstDeriv.y();

  if (denominator < EPSILON) {
    return {0};
  }

  double denom_double = CGAL::to_double(denominator);
  double power_result = std::pow(denom_double, 1.5);
  double num_double   = CGAL::to_double(numerator);
  return {num_double / power_result};
}

auto
NURBSCurve::torsion(Parameter parameter) const -> FT
{
  if (!is3D()) {
    return {0};
  }

  Point firstDeriv  = derivative(parameter, 1);
  Point secondDeriv = derivative(parameter, 2);
  Point thirdDeriv  = derivative(parameter, 3);

  FT crossX =
      firstDeriv.y() * secondDeriv.z() - firstDeriv.z() * secondDeriv.y();
  FT crossY =
      firstDeriv.z() * secondDeriv.x() - firstDeriv.x() * secondDeriv.z();
  FT crossZ =
      firstDeriv.x() * secondDeriv.y() - firstDeriv.y() * secondDeriv.x();

  FT scalarTripleProduct = crossX * thirdDeriv.x() + crossY * thirdDeriv.y() +
                           crossZ * thirdDeriv.z();
  FT crossMagnitudeSquared =
      crossX * crossX + crossY * crossY + crossZ * crossZ;

  if (crossMagnitudeSquared < EPSILON) {
    return {0};
  }

  return scalarTripleProduct / crossMagnitudeSquared;
}

auto
NURBSCurve::frenetFrame(Parameter parameter) const
    -> std::tuple<Point, Point, Point>
{
  Point tangentVec = tangent(parameter);

  if (is3D()) {
    Point binormalVec = binormal(parameter);
    Point normalVec   = Point(
        tangentVec.y() * binormalVec.z() - tangentVec.z() * binormalVec.y(),
        tangentVec.z() * binormalVec.x() - tangentVec.x() * binormalVec.z(),
        tangentVec.x() * binormalVec.y() - tangentVec.y() * binormalVec.x());
    return std::make_tuple(tangentVec, normalVec, binormalVec);
  }
  Point normalVec   = normal(parameter);
  Point binormalVec = Point(0, 0, 1);
  return std::make_tuple(tangentVec, normalVec, binormalVec);
}

auto
NURBSCurve::toLineString(unsigned int numSegments) const
    -> std::unique_ptr<LineString>
{
  auto lineString = std::make_unique<LineString>();
  if (_controlPoints.empty()) {
    return lineString;
  }

  auto bounds     = parameterBounds();
  FT   paramRange = bounds.second - bounds.first;

  for (unsigned int segIdx = 0; segIdx <= numSegments; ++segIdx) {
    FT parameter = bounds.first + (FT(segIdx) / FT(numSegments)) * paramRange;
    lineString->addPoint(evaluate(parameter));
  }

  return lineString;
}

auto
NURBSCurve::toLineStringAdaptive(FT tolerance, unsigned int minSegments,
                                 unsigned int maxSegments) const
    -> std::unique_ptr<LineString>
{
  auto lineString = std::make_unique<LineString>();
  if (_controlPoints.empty()) {
    return lineString;
  }

  auto bounds = parameterBounds();

  std::vector<Parameter> parameters;
  parameters.reserve(maxSegments + 1);

  FT deltaParam = (bounds.second - bounds.first) / FT(minSegments);
  for (unsigned int segIdx = 0; segIdx <= minSegments; ++segIdx) {
    parameters.push_back(bounds.first + FT(segIdx) * deltaParam);
  }

  bool refined = true;
  while (refined && parameters.size() <= maxSegments) {
    refined = false;

    for (size_t paramIdx = 0; paramIdx < parameters.size() - 1; ++paramIdx) {
      Parameter param1   = parameters[paramIdx];
      Parameter param2   = parameters[paramIdx + 1];
      Parameter midParam = (param1 + param2) / FT(2);

      Point point1   = evaluate(param1);
      Point point2   = evaluate(param2);
      Point midPoint = evaluate(midParam);

      Point interpolated = Point(
          (point1.x() + point2.x()) / FT(2), (point1.y() + point2.y()) / FT(2),
          point1.is3D() ? (point1.z() + point2.z()) / FT(2) : FT(0));

      FT deviation = algorithm::distance(midPoint, interpolated);

      if (deviation > tolerance) {
        parameters.insert(
            parameters.begin() +
                static_cast<std::vector<Parameter>::difference_type>(paramIdx +
                                                                     1),
            midParam);
        refined = true;
        break;
      }
    }
  }

  for (const auto &param : parameters) {
    lineString->addPoint(evaluate(param));
  }

  return lineString;
}

auto
NURBSCurve::parameterBounds() const -> std::pair<Parameter, Parameter>
{
  if (_knotVector.empty()) {
    return std::make_pair(FT(0), FT(1));
  }

  return std::make_pair(_knotVector[_degree],
                        _knotVector[_controlPoints.size()]);
}

auto
NURBSCurve::isPeriodic() const -> bool
{
  if (_controlPoints.size() < 2 || _degree == 0) {
    return false;
  }

  size_t minOverlap = _degree;
  if (_controlPoints.size() < 2 * minOverlap) {
    return false;
  }

  for (size_t idx = 0; idx < minOverlap; ++idx) {
    if (algorithm::distance(
            _controlPoints[idx],
            _controlPoints[_controlPoints.size() - minOverlap + idx]) >
        EPSILON) {
      return false;
    }
  }

  return true;
}

// NOLINTBEGIN(readability-function-cognitive-complexity)
auto
NURBSCurve::isPlanar(std::vector<FT> *plane) const -> bool
{
  if (!is3D() || _controlPoints.size() < 4) {
    if (plane != nullptr) {
      *plane = {FT(0), FT(0), FT(1), FT(0)};
    }
    return true;
  }

  for (size_t firstIdx = 0; firstIdx < _controlPoints.size() - 2; ++firstIdx) {
    for (size_t secondIdx = firstIdx + 1; secondIdx < _controlPoints.size() - 1;
         ++secondIdx) {
      for (size_t thirdIdx = secondIdx + 1; thirdIdx < _controlPoints.size();
           ++thirdIdx) {

        Point vec1 =
            Point(_controlPoints[secondIdx].x() - _controlPoints[firstIdx].x(),
                  _controlPoints[secondIdx].y() - _controlPoints[firstIdx].y(),
                  _controlPoints[secondIdx].z() - _controlPoints[firstIdx].z());

        Point vec2 =
            Point(_controlPoints[thirdIdx].x() - _controlPoints[firstIdx].x(),
                  _controlPoints[thirdIdx].y() - _controlPoints[firstIdx].y(),
                  _controlPoints[thirdIdx].z() - _controlPoints[firstIdx].z());

        Point normal = Point(vec1.y() * vec2.z() - vec1.z() * vec2.y(),
                             vec1.z() * vec2.x() - vec1.x() * vec2.z(),
                             vec1.x() * vec2.y() - vec1.y() * vec2.x());

        FT normalMagnitude = std::sqrt(
            CGAL::to_double(normal.x() * normal.x() + normal.y() * normal.y() +
                            normal.z() * normal.z()));

        if (normalMagnitude > EPSILON) {
          normal =
              Point(normal.x() / normalMagnitude, normal.y() / normalMagnitude,
                    normal.z() / normalMagnitude);

          FT planeD = -(normal.x() * _controlPoints[firstIdx].x() +
                        normal.y() * _controlPoints[firstIdx].y() +
                        normal.z() * _controlPoints[firstIdx].z());

          bool allPointsOnPlane = true;
          for (const auto &point : _controlPoints) {
            FT distance =
                CGAL::abs(normal.x() * point.x() + normal.y() * point.y() +
                          normal.z() * point.z() + planeD);
            if (distance > EPSILON) {
              allPointsOnPlane = false;
              break;
            }
          }

          if (allPointsOnPlane) {
            if (plane != nullptr) {
              *plane = {normal.x(), normal.y(), normal.z(), planeD};
            }
            return true;
          }
        }
      }
    }
  }

  return false;
}
// NOLINTEND(readability-function-cognitive-complexity)

auto
NURBSCurve::isLinear() const -> bool
{
  if (_controlPoints.size() < 2) {
    return true;
  }

  if (_degree > 1) {
    return false;
  }

  Point direction = Point(_controlPoints[1].x() - _controlPoints[0].x(),
                          _controlPoints[1].y() - _controlPoints[0].y(),
                          _controlPoints[1].is3D()
                              ? _controlPoints[1].z() - _controlPoints[0].z()
                              : FT(0));

  FT dirMagnitude = std::sqrt(CGAL::to_double(
      direction.x() * direction.x() + direction.y() * direction.y() +
      (direction.is3D() ? direction.z() * direction.z() : FT(0))));

  if (dirMagnitude < EPSILON) {
    return true;
  }

  direction = Point(direction.x() / dirMagnitude, direction.y() / dirMagnitude,
                    direction.is3D() ? direction.z() / dirMagnitude : FT(0));

  for (size_t ptIdx = 2; ptIdx < _controlPoints.size(); ++ptIdx) {
    Point vecToPoint =
        Point(_controlPoints[ptIdx].x() - _controlPoints[0].x(),
              _controlPoints[ptIdx].y() - _controlPoints[0].y(),
              _controlPoints[ptIdx].is3D()
                  ? _controlPoints[ptIdx].z() - _controlPoints[0].z()
                  : FT(0));

    FT vecMagnitude = std::sqrt(CGAL::to_double(
        vecToPoint.x() * vecToPoint.x() + vecToPoint.y() * vecToPoint.y() +
        (vecToPoint.is3D() ? vecToPoint.z() * vecToPoint.z() : FT(0))));

    if (vecMagnitude > EPSILON) {
      vecToPoint =
          Point(vecToPoint.x() / vecMagnitude, vecToPoint.y() / vecMagnitude,
                vecToPoint.is3D() ? vecToPoint.z() / vecMagnitude : FT(0));

      Point crossProduct = Point(
          direction.y() * vecToPoint.z() - direction.z() * vecToPoint.y(),
          direction.z() * vecToPoint.x() - direction.x() * vecToPoint.z(),
          direction.x() * vecToPoint.y() - direction.y() * vecToPoint.x());

      FT crossMagnitude =
          std::sqrt(CGAL::to_double(crossProduct.x() * crossProduct.x() +
                                    crossProduct.y() * crossProduct.y() +
                                    crossProduct.z() * crossProduct.z()));

      if (crossMagnitude > EPSILON) {
        return false;
      }
    }
  }

  return true;
}

auto
NURBSCurve::length(Parameter from, Parameter to, FT tolerance) const -> FT
{
  if (isEmpty()) {
    BOOST_THROW_EXCEPTION(Exception("Cannot compute length of empty curve"));
  }

  auto bounds = parameterBounds();

  if (from < FT(0)) {
    from = bounds.first;
  }
  if (to < FT(0)) {
    to = bounds.second;
  }

  if (from >= to) {
    return {0};
  }

  tolerance = CGAL::abs(tolerance);
  // NOLINTBEGIN(clang-analyzer-cplusplus.NewDeleteLeaks)
  if (tolerance <= FT(0)) {
    tolerance = EPSILON;
  }

  return computeArcLength(from, to, tolerance);
  // NOLINTEND(clang-analyzer-cplusplus.NewDeleteLeaks)
}

auto
NURBSCurve::parameterAtLength(FT arcLength, FT tolerance) const -> Parameter
{
  if (arcLength <= FT(0)) {
    return parameterBounds().first;
  }

  FT totalLength = length();
  if (arcLength >= totalLength) {
    return parameterBounds().second;
  }

  return findParameterByArcLength(arcLength, tolerance);
}

auto
NURBSCurve::reparameterizeByArcLength() const -> std::unique_ptr<Curve>
{
  if (isEmpty()) {
    return std::make_unique<NURBSCurve>();
  }

  FT totalLength = length();

  if (totalLength <= EPSILON) {
    return std::unique_ptr<NURBSCurve>(clone());
  }

  // Use a reasonable number of samples for reparameterization
  size_t numSamples =
      std::min(static_cast<size_t>(32),
               std::max(static_cast<size_t>(8), _controlPoints.size()));

  std::vector<Point> newControlPoints;
  newControlPoints.reserve(numSamples);

  for (size_t idx = 0; idx < numSamples; ++idx) {
    FT        targetLength = (FT(idx) / FT(numSamples - 1)) * totalLength;
    Parameter param        = parameterAtLength(targetLength, EPSILON);
    newControlPoints.push_back(evaluate(param));
  }

  // Use a reasonable degree for the reparameterized curve
  unsigned int newDegree =
      std::min(3U, static_cast<unsigned int>(newControlPoints.size() - 1));

  return std::make_unique<NURBSCurve>(newControlPoints, newDegree,
                                      KnotMethod::UNIFORM);
}

auto
NURBSCurve::split(Parameter parameter) const
    -> std::pair<std::unique_ptr<Curve>, std::unique_ptr<Curve>>
{
  if (isEmpty()) {
    return std::make_pair(std::make_unique<NURBSCurve>(),
                          std::make_unique<NURBSCurve>());
  }

  auto bounds = parameterBounds();

  parameter = std::max(bounds.first, std::min(bounds.second, parameter));

  if (parameter <= bounds.first) {
    auto emptyCurve = std::make_unique<NURBSCurve>();
    auto fullCurve  = std::unique_ptr<NURBSCurve>(clone());
    return std::make_pair(std::move(emptyCurve), std::move(fullCurve));
  }

  if (parameter >= bounds.second) {
    auto fullCurve  = std::unique_ptr<NURBSCurve>(clone());
    auto emptyCurve = std::make_unique<NURBSCurve>();
    return std::make_pair(std::move(fullCurve), std::move(emptyCurve));
  }

  auto firstCurve  = subcurve(bounds.first, parameter);
  auto secondCurve = subcurve(parameter, bounds.second);

  auto firstNurbs = std::unique_ptr<NURBSCurve>(
      dynamic_cast<NURBSCurve *>(firstCurve.release()));
  auto secondNurbs = std::unique_ptr<NURBSCurve>(
      dynamic_cast<NURBSCurve *>(secondCurve.release()));

  if (!firstNurbs) {
    firstNurbs = std::make_unique<NURBSCurve>();
  }
  if (!secondNurbs) {
    secondNurbs = std::make_unique<NURBSCurve>();
  }

  return std::make_pair(std::move(firstNurbs), std::move(secondNurbs));
}

auto
NURBSCurve::subcurve(Parameter from, Parameter to) const
    -> std::unique_ptr<Curve>
{
  if (isEmpty()) {
    return std::make_unique<NURBSCurve>();
  }

  if (from >= to) {
    return std::make_unique<NURBSCurve>();
  }

  auto bounds = parameterBounds();
  if (from < bounds.first) {
    from = bounds.first;
  }
  if (to > bounds.second) {
    to = bounds.second;
  }

  if (from >= to) {
    return std::make_unique<NURBSCurve>();
  }

  unsigned int numSamples =
      std::max(16U, static_cast<unsigned int>(_controlPoints.size()));

  std::vector<Point> subPoints;
  subPoints.reserve(numSamples + 1);

  FT paramRange = to - from;

  for (unsigned int sampleIdx = 0; sampleIdx <= numSamples; ++sampleIdx) {
    Parameter param = from + (FT(sampleIdx) / FT(numSamples)) * paramRange;
    Point     evaluatedPoint = evaluate(param);
    subPoints.push_back(evaluatedPoint);
  }

  unsigned int subDegree =
      std::min(_degree, static_cast<unsigned int>(subPoints.size() - 1));

  auto subKnots =
      generateKnotVector(subPoints, subDegree, KnotMethod::CHORD_LENGTH);

  std::vector<FT> subWeights;
  if (isRational()) {
    subWeights.resize(subPoints.size(), FT(1));
  }

  if (subWeights.empty()) {
    return std::make_unique<NURBSCurve>(subPoints, subDegree,
                                        KnotMethod::CHORD_LENGTH);
  }
  return std::make_unique<NURBSCurve>(subPoints, subWeights, subDegree,
                                      subKnots);
}

auto
NURBSCurve::reverse() const -> std::unique_ptr<Curve>
{
  if (_controlPoints.empty()) {
    return std::make_unique<NURBSCurve>();
  }

  std::vector<Point> reversedPoints(_controlPoints.rbegin(),
                                    _controlPoints.rend());
  std::vector<FT>    reversedWeights;

  if (!_weights.empty()) {
    reversedWeights.assign(_weights.rbegin(), _weights.rend());
  }

  std::vector<Knot> reversedKnots;
  if (!_knotVector.empty()) {
    reversedKnots.reserve(_knotVector.size());
    FT minKnot = _knotVector.front();
    FT maxKnot = _knotVector.back();

    // Proper knot reflection: new_knot = maxKnot + minKnot - old_knot
    for (auto knotIter = _knotVector.rbegin(); knotIter != _knotVector.rend();
         ++knotIter) {
      reversedKnots.push_back(maxKnot + minKnot - *knotIter);
    }
  }

  return std::make_unique<NURBSCurve>(reversedPoints, reversedWeights, _degree,
                                      reversedKnots);
}

auto
NURBSCurve::join(const Curve &other, Continuity /*continuity*/,
                 FT           tolerance) const -> std::unique_ptr<Curve>
{
  const auto *otherNurbs = dynamic_cast<const NURBSCurve *>(&other);
  if (otherNurbs == nullptr) {
    BOOST_THROW_EXCEPTION(
        Exception("Can only join NURBS curves with other NURBS curves"));
  }

  if (isEmpty() && otherNurbs->isEmpty()) {
    return std::make_unique<NURBSCurve>();
  }

  if (isEmpty()) {
    return std::unique_ptr<NURBSCurve>(otherNurbs->clone());
  }

  if (otherNurbs->isEmpty()) {
    return std::unique_ptr<NURBSCurve>(clone());
  }

  Point thisEnd    = endPoint();
  Point otherStart = otherNurbs->startPoint();

  if (algorithm::distance(thisEnd, otherStart) > tolerance) {
    BOOST_THROW_EXCEPTION(
        Exception("Curves are not adjacent within tolerance"));
  }

  // Prepare joined control points
  std::vector<Point> joinedPoints = _controlPoints;
  const auto        &otherPoints  = otherNurbs->controlPoints();

  // Skip the first control point of the second curve to avoid duplication
  joinedPoints.insert(joinedPoints.end(), otherPoints.begin() + 1,
                      otherPoints.end());

  // Prepare joined weights
  std::vector<FT> joinedWeights;
  if (!_weights.empty() || !otherNurbs->weights().empty()) {
    // If either curve has weights, we need to handle weights for both
    const std::vector<FT> &firstWeights =
        _weights.empty() ? std::vector<FT>(_controlPoints.size(), FT(1))
                         : _weights;
    const std::vector<FT> &secondWeights =
        otherNurbs->weights().empty()
            ? std::vector<FT>(otherNurbs->numControlPoints(), FT(1))
            : otherNurbs->weights();

    joinedWeights = firstWeights;
    // Skip the first weight of the second curve (corresponding to the first
    // control point)
    joinedWeights.insert(joinedWeights.end(), secondWeights.begin() + 1,
                         secondWeights.end());
  }

  // Use the higher degree for the joined curve
  unsigned int joinedDegree = std::max(_degree, otherNurbs->degree());
  joinedDegree              = std::min(joinedDegree,
                                       static_cast<unsigned int>(joinedPoints.size() - 1));

  // Properly merge knot vectors
  std::vector<FT> joinedKnots;
  if (!_knotVector.empty() && !otherNurbs->knotVector().empty()) {
    // Get the knot vectors
    const auto &firstKnots  = _knotVector;
    const auto &secondKnots = otherNurbs->knotVector();

    // Add the first curve's knots
    joinedKnots = firstKnots;

    // Find the last parameter value of the first curve
    FT lastParam = firstKnots.back();

    // Add the second curve's knots, shifted and with first multiplicity removed
    for (size_t i = joinedDegree + 1; i < secondKnots.size(); ++i) {
      joinedKnots.push_back(lastParam + secondKnots[i] - secondKnots[0]);
    }
  } else {
    // Generate appropriate knot vector for the joined curve
    joinedKnots =
        generateKnotVector(joinedPoints, joinedDegree, KnotMethod::UNIFORM);
  }

  if (joinedWeights.empty()) {
    // Use uniform weights for non-rational case
    std::vector<FT> uniformWeights(joinedPoints.size(), FT(1));
    return std::make_unique<NURBSCurve>(joinedPoints, uniformWeights,
                                        joinedDegree, joinedKnots);
  }
  return std::make_unique<NURBSCurve>(joinedPoints, joinedWeights, joinedDegree,
                                      joinedKnots);
}

auto
NURBSCurve::closestPoint(const Point &point, Parameter *outParameter) const
    -> Point
{
  if (_controlPoints.empty()) {
    BOOST_THROW_EXCEPTION(
        Exception("Cannot find closest point on empty curve"));
  }

  auto [closestPt, parameter] = projectPointToCurve(point);

  if (outParameter != nullptr) {
    *outParameter = parameter;
  }

  return closestPt;
}

void
NURBSCurve::setControlPoints(const std::vector<Point> &controlPoints)
{
  _controlPoints = controlPoints;

  if (!_weights.empty() && _weights.size() != controlPoints.size()) {
    _weights.resize(controlPoints.size(), FT(1));
  }

  // Only generate knot vector when none exists and we have control points
  if (_knotVector.empty() && !controlPoints.empty()) {
    _knotVector =
        generateKnotVector(controlPoints, _degree, KnotMethod::CHORD_LENGTH);
  }

  auto [isValid, reason] = validateData();
  if (!isValid) {
    BOOST_THROW_EXCEPTION(Exception("Invalid control points: " + reason));
  }
}

auto
NURBSCurve::controlPointN(size_t index) -> Point &
{
  if (index >= _controlPoints.size()) {
    BOOST_THROW_EXCEPTION(Exception("Control point index out of bounds"));
  }
  return _controlPoints[index];
}

auto
NURBSCurve::controlPointN(size_t index) const -> const Point &
{
  if (index >= _controlPoints.size()) {
    BOOST_THROW_EXCEPTION(Exception("Control point index out of bounds"));
  }
  return _controlPoints[index];
}

void
NURBSCurve::setControlPoint(size_t index, const Point &point)
{
  if (index >= _controlPoints.size()) {
    BOOST_THROW_EXCEPTION(Exception("Control point index out of bounds"));
  }

  if (!_controlPoints.empty() &&
      (_controlPoints[0].is3D() != point.is3D() ||
       _controlPoints[0].isMeasured() != point.isMeasured())) {
    BOOST_THROW_EXCEPTION(Exception("Control point dimension mismatch"));
  }

  _controlPoints[index] = point;
}

void
NURBSCurve::setWeights(const std::vector<FT> &weights)
{
  if (!weights.empty() && weights.size() != _controlPoints.size()) {
    BOOST_THROW_EXCEPTION(
        Exception("Weight count (" + std::to_string(weights.size()) +
                  ") does not match control point count (" +
                  std::to_string(_controlPoints.size()) + ")"));
  }

  for (size_t weightIdx = 0; weightIdx < weights.size(); ++weightIdx) {
    if (weights[weightIdx] <= FT(0)) {
      BOOST_THROW_EXCEPTION(Exception("Weight at index " +
                                      std::to_string(weightIdx) +
                                      " must be positive"));
    }
  }

  _weights = weights;
}

auto
NURBSCurve::weight(size_t index) const -> FT
{
  if (_weights.empty()) {
    return {1};
  }
  if (index >= _weights.size()) {
    BOOST_THROW_EXCEPTION(Exception("Weight index out of bounds"));
  }
  return _weights[index];
}

void
NURBSCurve::setWeight(size_t index, const FT &weight)
{
  if (weight <= FT(0)) {
    BOOST_THROW_EXCEPTION(Exception("Weight must be positive"));
  }

  if (_weights.empty()) {
    _weights.resize(_controlPoints.size(), FT(1));
  }

  if (index >= _weights.size()) {
    BOOST_THROW_EXCEPTION(Exception("Weight index out of bounds"));
  }
  _weights[index] = weight;
}

auto
NURBSCurve::isRational() const -> bool
{
  if (_weights.empty() || _weights.size() < 2) {
    return false;
  }

  FT       firstWeight = _weights[0];
  const FT tolerance   = EPSILON;

  return std::any_of(_weights.begin(), _weights.end(),
                     [firstWeight, tolerance](const auto &weight) {
                       return CGAL::abs(weight - firstWeight) > tolerance;
                     });
}

auto
NURBSCurve::isBezier() const -> bool
{
  if (_knotVector.empty() || _controlPoints.empty()) {
    return false;
  }

  size_t expectedSize = static_cast<size_t>(2) * (_degree + 1);
  if (_knotVector.size() != expectedSize) {
    return false;
  }

  for (unsigned int idx = 0; idx <= _degree; ++idx) {
    if (CGAL::abs(_knotVector[idx] - _knotVector[0]) > EPSILON) {
      return false;
    }
  }

  for (unsigned int idx = _degree + 1; idx < _knotVector.size(); ++idx) {
    if (CGAL::abs(_knotVector[idx] - _knotVector.back()) > EPSILON) {
      return false;
    }
  }

  return true;
}

auto
NURBSCurve::isBSpline() const -> bool
{
  return !isRational();
}

void
NURBSCurve::setKnotVector(const std::vector<Knot> &knots)
{
  _knotVector            = knots;
  auto [isValid, reason] = validateData();
  if (!isValid) {
    BOOST_THROW_EXCEPTION(Exception("Invalid knot vector: " + reason));
  }
}

auto
NURBSCurve::knot(size_t index) const -> Knot
{
  if (index >= _knotVector.size()) {
    BOOST_THROW_EXCEPTION(Exception("Knot index out of bounds"));
  }
  return _knotVector[index];
}

auto
NURBSCurve::knotMultiplicity(const Knot &value, const FT &tolerance) const
    -> unsigned int
{
  unsigned int multiplicity = 0;
  for (const auto &knot : _knotVector) {
    if (CGAL::abs(knot - value) <= tolerance) {
      ++multiplicity;
    }
  }
  return multiplicity;
}

// NOLINTBEGIN(readability-convert-member-functions-to-static)
auto
NURBSCurve::insertKnot(const Knot &parameter, unsigned int times) const
    -> std::unique_ptr<NURBSCurve>
{
  BOOST_THROW_EXCEPTION(
      Exception("insertKnot: NURBS knot insertion not yet implemented"));
}

auto
NURBSCurve::refineKnotVector(const std::vector<Knot> &newKnots) const
    -> std::unique_ptr<NURBSCurve>
{
  BOOST_THROW_EXCEPTION(Exception(
      "refineKnotVector: NURBS knot vector refinement not yet implemented"));
}

auto
NURBSCurve::elevateDegree(unsigned int times) const
    -> std::unique_ptr<NURBSCurve>
{
  BOOST_THROW_EXCEPTION(
      Exception("elevateDegree: NURBS degree elevation not yet implemented"));
}

auto
NURBSCurve::reduceDegree(const FT &tolerance) const
    -> std::unique_ptr<NURBSCurve>
{
  BOOST_THROW_EXCEPTION(
      Exception("reduceDegree: NURBS degree reduction not yet implemented"));
}

auto
NURBSCurve::removeKnots(const FT &tolerance) const
    -> std::unique_ptr<NURBSCurve>
{
  BOOST_THROW_EXCEPTION(
      Exception("removeKnots: NURBS knot removal not yet implemented"));
}
// NOLINTEND(readability-convert-member-functions-to-static)

auto
NURBSCurve::validateData() const -> std::pair<bool, std::string>
{
  return validateNURBSData(_controlPoints, _weights, _degree, _knotVector);
}

// NOLINTBEGIN(readability-function-cognitive-complexity)
auto
NURBSCurve::validateNURBSData(const std::vector<Point> &controlPoints,
                              const std::vector<FT>    &weights,
                              unsigned int              degree,
                              const std::vector<Knot>  &knots)
    -> std::pair<bool, std::string>
{
  if (controlPoints.empty()) {
    if (!knots.empty() || !weights.empty()) {
      return {false, "Empty curve with non-empty knot vector or weights"};
    }
    return {true, ""};
  }

  for (const auto &point : controlPoints) {
    if ((controlPoints[0].is3D() != point.is3D() ||
         controlPoints[0].isMeasured() != point.isMeasured())) {
      return {false, "Control point dimension mismatch"};
    }
  }

  if (degree >= controlPoints.size()) {
    return {false, "Degree must be less than control points count"};
  }

  const size_t expectedKnots = controlPoints.size() + degree + 1;
  if (knots.size() != expectedKnots) {
    return {false, "Invalid knot vector size"};
  }

  for (size_t idx = 1; idx < knots.size(); ++idx) {
    if (knots[idx] < knots[idx - 1]) {
      return {false, "Knot vector is not non-decreasing"};
    }
  }

  if (!weights.empty()) {
    if (weights.size() != controlPoints.size()) {
      return {false, "Weight count does not match control point count"};
    }

    for (const auto &weight : weights) {
      if (weight <= FT(0) || !CGAL::is_finite(weight)) {
        return {false, "Invalid weight value"};
      }
    }
  }

  if (!knots.empty() && !controlPoints.empty()) {
    const FT &startParam = knots[degree];
    const FT &endParam   = knots[controlPoints.size()];

    if (startParam > endParam) {
      return {false, "Invalid parameter bounds"};
    }
  }

  return {true, ""};
}
// NOLINTEND(readability-function-cognitive-complexity)

auto
NURBSCurve::getCurveStatistics() const -> std::map<std::string, double>
{
  std::map<std::string, double> stats;

  stats["num_control_points"] = static_cast<double>(_controlPoints.size());
  stats["degree"]             = static_cast<double>(_degree);
  stats["num_knots"]          = static_cast<double>(_knotVector.size());
  stats["is_rational"]        = isRational() ? 1.0 : 0.0;
  stats["is_bezier"]          = isBezier() ? 1.0 : 0.0;
  stats["is_3d"]              = is3D() ? 1.0 : 0.0;
  stats["is_measured"]        = isMeasured() ? 1.0 : 0.0;

  if (!isEmpty()) {
    stats["curve_length"] = CGAL::to_double(length());

    auto bounds        = parameterBounds();
    stats["param_min"] = CGAL::to_double(bounds.first);
    stats["param_max"] = CGAL::to_double(bounds.second);
  }

  return stats;
}

auto
NURBSCurve::findSpan(const Parameter &parameter) const -> size_t
{
  // Delegate to the shared implementation
  return findKnotSpanImpl(parameter, _degree, _knotVector,
                          _controlPoints.size());
}

// NOLINTBEGIN(readability-function-cognitive-complexity)
auto
NURBSCurve::deBoorRational(size_t span, const Parameter &parameter) const
    -> Point
{
  // De Boor's algorithm for rational B-spline evaluation
  // Modified to handle weights and XYZM coordinates properly
  // Based on Algorithm A4.1 from "The NURBS Book"
  struct HomogeneousPoint {
    FT     weightedX, weightedY, weightedZ, weight;
    double measure;

    HomogeneousPoint()
        : weightedX(0), weightedY(0), weightedZ(0), weight(1), measure(0)
    {
    }

    HomogeneousPoint(const Point &point, const FT &pointWeight)
        : weightedX(point.x() * pointWeight),
          weightedY(point.y() * pointWeight),
          weightedZ(point.is3D() ? point.z() * pointWeight : FT(0)),
          weight(pointWeight), measure(point.isMeasured() ? point.m() : 0.0)
    {
    }
  };

  std::vector<HomogeneousPoint> temp(_degree + 1);

  for (unsigned int controlIdx = 0; controlIdx <= _degree; ++controlIdx) {
    size_t cpIndex = span - _degree + controlIdx;
    if (cpIndex >= _controlPoints.size()) {
      cpIndex = _controlPoints.size() - 1;
    }
    FT pointWeight   = _weights.empty() ? FT(1) : _weights[cpIndex];
    temp[controlIdx] = HomogeneousPoint(_controlPoints[cpIndex], pointWeight);
  }

  for (unsigned int level = 1; level <= _degree; ++level) {
    for (unsigned int basisIdx = _degree; basisIdx >= level; --basisIdx) {
      size_t knotIdx = span - _degree + basisIdx;

      if (knotIdx + _degree - level + 1 >= _knotVector.size() ||
          knotIdx >= _knotVector.size()) {
        continue;
      }

      FT denominator =
          _knotVector[knotIdx + _degree - level + 1] - _knotVector[knotIdx];
      FT alpha = FT(0);

      if (CGAL::abs(denominator) > EPSILON) {
        alpha = (parameter - _knotVector[knotIdx]) / denominator;
      } else {
        // When denominator is zero (repeated knots), use proper limit
        // In this case, alpha should be 0 if parameter <
        // knotVector[knotIdx+degree-level+1] and 1 if parameter >=
        // knotVector[knotIdx+degree-level+1]
        if (parameter >= _knotVector[knotIdx + _degree - level + 1]) {
          alpha = FT(1);
        } else {
          alpha = FT(0);
        }
      }

      alpha = std::max(FT(0), std::min(FT(1), alpha));

      temp[basisIdx].weightedX =
          (FT(1) - alpha) * temp[basisIdx - 1].weightedX +
          alpha * temp[basisIdx].weightedX;
      temp[basisIdx].weightedY =
          (FT(1) - alpha) * temp[basisIdx - 1].weightedY +
          alpha * temp[basisIdx].weightedY;
      temp[basisIdx].weightedZ =
          (FT(1) - alpha) * temp[basisIdx - 1].weightedZ +
          alpha * temp[basisIdx].weightedZ;
      temp[basisIdx].weight = (FT(1) - alpha) * temp[basisIdx - 1].weight +
                              alpha * temp[basisIdx].weight;

      if (isMeasured()) {
        double alphaDouble = CGAL::to_double(alpha);
        temp[basisIdx].measure =
            (1.0 - alphaDouble) * temp[basisIdx - 1].measure +
            alphaDouble * temp[basisIdx].measure;
      }
    }
  }

  HomogeneousPoint &result = temp[_degree];

  if (CGAL::abs(result.weight) < EPSILON) {
    BOOST_THROW_EXCEPTION(Exception("Division by zero in NURBS evaluation"));
  }

  FT invWeight = FT(1) / result.weight;

  if (is3D() && isMeasured()) {
    return {result.weightedX * invWeight, result.weightedY * invWeight,
            result.weightedZ * invWeight, result.measure};
  }
  if (is3D()) {
    return {result.weightedX * invWeight, result.weightedY * invWeight,
            result.weightedZ * invWeight};
  }
  if (isMeasured()) {
    Point point{result.weightedX * invWeight, result.weightedY * invWeight};
    point.setM(result.measure);
    return point;
  }
  return {result.weightedX * invWeight, result.weightedY * invWeight};
}
// NOLINTEND(readability-function-cognitive-complexity)

// NOLINTBEGIN(readability-function-cognitive-complexity)
auto
NURBSCurve::computeDerivatives(const Parameter &parameter,
                               unsigned int     maxOrder) const
    -> std::vector<Point>
{
  std::vector<Point> derivatives;
  derivatives.reserve(maxOrder + 1);

  derivatives.push_back(evaluate(parameter));

  if (maxOrder > 0) {
    if (_degree == 1 && _controlPoints.size() == 2) {
      Point derivative =
          Point(_controlPoints[1].x() - _controlPoints[0].x(),
                _controlPoints[1].y() - _controlPoints[0].y(),
                _controlPoints[1].is3D()
                    ? _controlPoints[1].z() - _controlPoints[0].z()
                    : FT(0));
      derivatives.push_back(derivative);

      for (unsigned int orderIdx = 2; orderIdx <= maxOrder; ++orderIdx) {
        if (is3D() && isMeasured()) {
          derivatives.emplace_back(FT(0), FT(0), FT(0), 0.0);
        } else if (is3D()) {
          derivatives.emplace_back(FT(0), FT(0), FT(0));
        } else if (isMeasured()) {
          Point zeroPoint(FT(0), FT(0));
          zeroPoint.setM(0.0);
          derivatives.push_back(zeroPoint);
        } else {
          derivatives.emplace_back(FT(0), FT(0));
        }
      }
    } else {
      auto bounds = parameterBounds();

      for (unsigned int orderIdx = 1; orderIdx <= maxOrder; ++orderIdx) {
        FT param1     = std::max(bounds.first, parameter - EPSILON);
        FT param2     = std::min(bounds.second, parameter + EPSILON);
        FT deltaParam = param2 - param1;

        if (deltaParam > FT(0)) {
          Point point1 = evaluate(param1);
          Point point2 = evaluate(param2);

          Point derivative(
              (point2.x() - point1.x()) / deltaParam,
              (point2.y() - point1.y()) / deltaParam,
              point2.is3D() ? (point2.z() - point1.z()) / deltaParam : FT(0));

          if (isMeasured()) {
            derivative.setM((point2.m() - point1.m()) /
                            CGAL::to_double(deltaParam));
          }

          derivatives.push_back(derivative);
        } else {
          if (is3D() && isMeasured()) {
            derivatives.emplace_back(FT(0), FT(0), FT(0), 0.0);
          } else if (is3D()) {
            derivatives.emplace_back(FT(0), FT(0), FT(0));
          } else if (isMeasured()) {
            Point zeroPoint(FT(0), FT(0));
            zeroPoint.setM(0.0);
            derivatives.push_back(zeroPoint);
          } else {
            derivatives.emplace_back(FT(0), FT(0));
          }
        }
      }
    }
  }

  return derivatives;
}
// NOLINTEND(readability-function-cognitive-complexity)

auto
NURBSCurve::generateKnotVector(const std::vector<Point> &points,
                               unsigned int degree, KnotMethod method)
    -> std::vector<Knot>
{
  std::vector<Knot> knots;
  size_t            numPoints = points.size();
  size_t            numKnots  = numPoints + degree + 1;

  if (numPoints == 0) {
    return knots;
  }

  if (degree == 0) {
    knots.reserve(numKnots);
    for (size_t pointIdx = 0; pointIdx < numPoints; ++pointIdx) {
      knots.emplace_back(pointIdx);
    }
    knots.emplace_back(numPoints);
    return knots;
  }

  knots.reserve(numKnots);

  if (method == KnotMethod::UNIFORM) {
    for (unsigned int degreeIdx = 0; degreeIdx <= degree; ++degreeIdx) {
      knots.emplace_back(0);
    }

    size_t numInternal = numKnots - (static_cast<size_t>(2) * (degree + 1));
    for (size_t internalIdx = 1; internalIdx <= numInternal; ++internalIdx) {
      knots.push_back(FT(internalIdx) / FT(numInternal + 1));
    }

    for (unsigned int degreeIdx = 0; degreeIdx <= degree; ++degreeIdx) {
      knots.emplace_back(1);
    }
  } else {
    auto parameters = computeParameters(points, method);

    for (unsigned int degreeIdx = 0; degreeIdx <= degree; ++degreeIdx) {
      knots.push_back(parameters.front());
    }

    if (numPoints > degree + 1) {
      for (size_t pointIdx = 1; pointIdx < numPoints - degree; ++pointIdx) {
        FT sum = FT(0);
        for (unsigned int degreeIdx = 0; degreeIdx < degree; ++degreeIdx) {
          sum += parameters[pointIdx + degreeIdx];
        }
        knots.push_back(sum / FT(degree));
      }
    }

    for (unsigned int degreeIdx = 0; degreeIdx <= degree; ++degreeIdx) {
      knots.push_back(parameters.back());
    }
  }

  return knots;
}

auto
NURBSCurve::computeParameters(const std::vector<Point> &points,
                              KnotMethod method) -> std::vector<Parameter>
{
  std::vector<Parameter> parameters;
  parameters.reserve(points.size());

  if (points.empty()) {
    return parameters;
  }

  if (method == KnotMethod::UNIFORM) {
    for (size_t pointIdx = 0; pointIdx < points.size(); ++pointIdx) {
      parameters.push_back(FT(pointIdx) / FT(points.size() - 1));
    }
  } else {
    parameters.emplace_back(0);
    FT totalLength = FT(0);

    std::vector<FT> distances;
    distances.reserve(points.size() - 1);

    for (size_t pointIdx = 1; pointIdx < points.size(); ++pointIdx) {
      FT distance = algorithm::distance(points[pointIdx - 1], points[pointIdx]);

      if (method == KnotMethod::CENTRIPETAL) {
        distance = FT(std::sqrt(CGAL::to_double(distance)));
      }

      distances.push_back(distance);
      totalLength += distance;
      parameters.push_back(totalLength);
    }

    if (totalLength > FT(0)) {
      for (auto &param : parameters) {
        param /= totalLength;
      }
    }
  }

  return parameters;
}

auto
NURBSCurve::checkDimensionalConsistency() const -> bool
{
  if (_controlPoints.empty()) {
    return true;
  }

  bool first3D       = _controlPoints[0].is3D();
  bool firstMeasured = _controlPoints[0].isMeasured();

  return std::all_of(_controlPoints.begin(), _controlPoints.end(),
                     [first3D, firstMeasured](const auto &point) {
                       return point.is3D() == first3D &&
                              point.isMeasured() == firstMeasured;
                     });
}

auto
NURBSCurve::computeArcLength(const Parameter &startParam,
                             const Parameter &endParam,
                             const FT        &tolerance) const -> FT
{
  if (startParam >= endParam) {
    return {0};
  }

  /**
   * Arc length computation using adaptive Simpson quadrature.
   *
   * Integrates the speed function ||C'(t)|| over [startParam, endParam]:
   * L = ∫[a,b] ||dC/dt|| dt
   *
   * The adaptive approach subdivides intervals where Simpson's rule
   * and its halving disagree by more than 15*tolerance, ensuring
   * accurate results for curves with varying curvature.
   */
  return adaptiveSimpsonArcLength(startParam, endParam, tolerance);
}

auto
NURBSCurve::speedFunction(Parameter t) const -> FT
{
  // Compute ||C'(t)|| - the magnitude of the derivative
  Point derivative = this->derivative(std::move(t), 1);

  if (is3D()) {
    FT dx = derivative.x();
    FT dy = derivative.y();
    FT dz = derivative.z();
    // Use CGAL::to_double for sqrt calculation, then convert back to FT
    double magnitude = std::sqrt(CGAL::to_double(dx * dx + dy * dy + dz * dz));
    return {magnitude};
  }
  FT dx = derivative.x();
  FT dy = derivative.y();
  // Use CGAL::to_double for sqrt calculation, then convert back to FT
  double magnitude = std::sqrt(CGAL::to_double(dx * dx + dy * dy));
  return {magnitude};
}

auto
NURBSCurve::simpsonRule(const Parameter &a, const Parameter &b) const -> FT
{
  // Basic Simpson's rule: (b-a)/6 * [f(a) + 4*f((a+b)/2) + f(b)]
  Parameter mid  = (a + b) / FT(2);
  FT        fa   = speedFunction(a);
  FT        fmid = speedFunction(mid);
  FT        fb   = speedFunction(b);

  return (b - a) / FT(6) * (fa + FT(4) * fmid + fb);
}

auto
NURBSCurve::adaptiveSimpsonArcLength(const Parameter &a, const Parameter &b,
                                     const FT    &tolerance,
                                     unsigned int maxDepth) const -> FT
{
  if (maxDepth == 0) {
    // Fallback to basic Simpson's rule when max recursion reached
    return simpsonRule(a, b);
  }

  Parameter mid = (a + b) / FT(2);

  // Compute Simpson's rule for whole interval and two halves
  FT wholeInterval = simpsonRule(a, b);
  FT leftHalf      = simpsonRule(a, mid);
  FT rightHalf     = simpsonRule(mid, b);
  FT twoHalves     = leftHalf + rightHalf;

  // Error estimate: |S(a,b) - S(a,m) - S(m,b)| / 15
  FT error = CGAL::abs(wholeInterval - twoHalves) / FT(15);

  if (error <= tolerance) {
    // Richardson extrapolation: more accurate result
    return twoHalves + (twoHalves - wholeInterval) / FT(15);
  } // Recursively subdivide with tighter tolerance
  FT halfTolerance = tolerance / FT(2);
  return adaptiveSimpsonArcLength(a, mid, halfTolerance, maxDepth - 1) +
         adaptiveSimpsonArcLength(mid, b, halfTolerance, maxDepth - 1);
}

auto
NURBSCurve::findParameterByArcLength(const FT &targetLength,
                                     const FT &tolerance) const -> Parameter
{
  auto bounds = parameterBounds();

  // Use binary search with fixed iterations to avoid infinite loops
  Parameter low  = bounds.first;
  Parameter high = bounds.second;

  const unsigned int maxIterations =
      50; // Limit iterations to prevent infinite loops

  for (unsigned int iter = 0; iter < maxIterations; ++iter) {
    Parameter mid         = (low + high) / FT(2);
    FT        lengthAtMid = computeArcLength(bounds.first, mid, tolerance);

    if (CGAL::abs(lengthAtMid - targetLength) <= tolerance) {
      return mid;
    }

    if (lengthAtMid < targetLength) {
      low = mid;
    } else {
      high = mid;
    }

    // Check for convergence
    if (high - low < EPSILON) {
      break;
    }
  }

  return (low + high) / FT(2);
}

auto
NURBSCurve::projectPointToCurve(const Point &point, const FT &tolerance,
                                unsigned int maxIterations) const
    -> std::pair<Point, Parameter>
{
  if (_controlPoints.empty()) {
    BOOST_THROW_EXCEPTION(Exception("Cannot project to empty curve"));
  }

  FT        minDistance = algorithm::distance(point, _controlPoints[0]);
  Parameter bestParam   = parameterBounds().first;

  auto bounds     = parameterBounds();
  FT   paramRange = bounds.second - bounds.first;

  unsigned int numSamples = 50;
  for (unsigned int sampleIdx = 0; sampleIdx <= numSamples; ++sampleIdx) {
    Parameter param =
        bounds.first + (FT(sampleIdx) / FT(numSamples)) * paramRange;
    Point curvePoint = evaluate(param);
    FT    distance   = algorithm::distance(point, curvePoint);

    if (distance < minDistance) {
      minDistance = distance;
      bestParam   = param;
    }
  }

  for (unsigned int iteration = 0; iteration < maxIterations; ++iteration) {
    Point curvePoint = evaluate(bestParam);
    Point firstDeriv = derivative(bestParam, 1);

    Point diff = Point(curvePoint.x() - point.x(), curvePoint.y() - point.y(),
                       curvePoint.is3D() ? curvePoint.z() - point.z() : FT(0));

    FT numerator = diff.x() * firstDeriv.x() + diff.y() * firstDeriv.y() +
                   (diff.is3D() ? diff.z() * firstDeriv.z() : FT(0));

    FT denominator =
        firstDeriv.x() * firstDeriv.x() + firstDeriv.y() * firstDeriv.y() +
        (firstDeriv.is3D() ? firstDeriv.z() * firstDeriv.z() : FT(0));

    if (denominator < EPSILON) {
      break;
    }

    FT paramUpdate = numerator / denominator;

    if (CGAL::abs(paramUpdate) <= tolerance) {
      break;
    }

    bestParam -= paramUpdate;
    bestParam = std::max(bounds.first, std::min(bounds.second, bestParam));
  }

  return std::make_pair(evaluate(bestParam), bestParam);
}

} // namespace SFCGAL
