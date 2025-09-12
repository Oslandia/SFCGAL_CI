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
using FT = NURBSCurve::FT;

//-- Constructors and initialization

NURBSCurve::NURBSCurve() : _degree(0), _fitTolerance(FT(0)) {}

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

//-- Factory methods

auto
NURBSCurve::fromBezier(const std::vector<Point> &controlPoints)
    -> std::unique_ptr<NURBSCurve>
{
  if (controlPoints.empty()) {
    return std::make_unique<NURBSCurve>();
  }

  unsigned int degree = static_cast<unsigned int>(controlPoints.size() - 1);

  std::vector<Knot> knots;
  knots.reserve(2 * (degree + 1));

  for (unsigned int idx = 0; idx <= degree; ++idx) {
    knots.push_back(FT(0));
  }
  for (unsigned int idx = 0; idx <= degree; ++idx) {
    knots.push_back(FT(1));
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

auto
NURBSCurve::createCircularArc(const Point &center, FT radius, FT startAngle,
                              FT endAngle, const Point &normal)
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
  if (angleSpan <= M_PI + 1e-10) {
    double halfSpan    = angleSpan / 2.0;
    double midAngle    = (startAngleDouble + endAngleDouble) / 2.0;
    double cosHalfSpan = std::cos(halfSpan);

    // Start point
    controlPoints.emplace_back(
        center.x() + FT(radiusDouble * std::cos(startAngleDouble)),
        center.y() + FT(radiusDouble * std::sin(startAngleDouble)),
        center.is3D() ? center.z() : FT(0),
        center.isMeasured() ? center.m() : NaN(), dimType);
    weights.push_back(FT(1.0));

    // Middle control point (off the circle)
    controlPoints.emplace_back(
        center.x() + FT(radiusDouble * std::cos(midAngle) / cosHalfSpan),
        center.y() + FT(radiusDouble * std::sin(midAngle) / cosHalfSpan),
        center.is3D() ? center.z() : FT(0),
        center.isMeasured() ? center.m() : NaN(), dimType);
    weights.push_back(FT(cosHalfSpan));

    // End point
    controlPoints.emplace_back(
        center.x() + FT(radiusDouble * std::cos(endAngleDouble)),
        center.y() + FT(radiusDouble * std::sin(endAngleDouble)),
        center.is3D() ? center.z() : FT(0),
        center.isMeasured() ? center.m() : NaN(), dimType);
    weights.push_back(FT(1.0));

    std::vector<Knot> knots = {FT(0), FT(0), FT(0), FT(1), FT(1), FT(1)};
    return std::make_unique<NURBSCurve>(controlPoints, weights, 2, knots);
  }

  // For arcs > π, split into two segments
  double midAngle = startAngleDouble + angleSpan / 2.0;

  // First segment: start to mid
  double halfSpan1    = (midAngle - startAngleDouble) / 2.0;
  double midAngle1    = (startAngleDouble + midAngle) / 2.0;
  double cosHalfSpan1 = std::cos(halfSpan1);

  controlPoints.emplace_back(
      center.x() + FT(radiusDouble * std::cos(startAngleDouble)),
      center.y() + FT(radiusDouble * std::sin(startAngleDouble)),
      center.is3D() ? center.z() : FT(0),
      center.isMeasured() ? center.m() : NaN(), dimType);
  weights.push_back(FT(1.0));

  controlPoints.emplace_back(
      center.x() + FT(radiusDouble * std::cos(midAngle1) / cosHalfSpan1),
      center.y() + FT(radiusDouble * std::sin(midAngle1) / cosHalfSpan1),
      center.is3D() ? center.z() : FT(0),
      center.isMeasured() ? center.m() : NaN(), dimType);
  weights.push_back(FT(cosHalfSpan1));

  controlPoints.emplace_back(center.x() + FT(radiusDouble * std::cos(midAngle)),
                             center.y() + FT(radiusDouble * std::sin(midAngle)),
                             center.is3D() ? center.z() : FT(0),
                             center.isMeasured() ? center.m() : NaN(), dimType);
  weights.push_back(FT(1.0));

  // Second segment: mid to end
  double halfSpan2    = (endAngleDouble - midAngle) / 2.0;
  double midAngle2    = (midAngle + endAngleDouble) / 2.0;
  double cosHalfSpan2 = std::cos(halfSpan2);

  controlPoints.emplace_back(
      center.x() + FT(radiusDouble * std::cos(midAngle2) / cosHalfSpan2),
      center.y() + FT(radiusDouble * std::sin(midAngle2) / cosHalfSpan2),
      center.is3D() ? center.z() : FT(0),
      center.isMeasured() ? center.m() : NaN(), dimType);
  weights.push_back(FT(cosHalfSpan2));

  controlPoints.emplace_back(
      center.x() + FT(radiusDouble * std::cos(endAngleDouble)),
      center.y() + FT(radiusDouble * std::sin(endAngleDouble)),
      center.is3D() ? center.z() : FT(0),
      center.isMeasured() ? center.m() : NaN(), dimType);
  weights.push_back(FT(1.0));

  // Knot vector for two Bézier segments
  std::vector<Knot> knots = {FT(0),   FT(0), FT(0), FT(0.5),
                             FT(0.5), FT(1), FT(1), FT(1)};

  return std::make_unique<NURBSCurve>(controlPoints, weights, 2, knots);
}

// Helper functions for different end condition interpolations

// Static helper functions for matrix-based interpolation
static auto
findKnotSpan(FT parameter, unsigned int degree, const std::vector<FT> &knots)
    -> size_t
{
  size_t n = knots.size() - degree - 1;

  if (parameter >= knots[n]) {
    return n - 1;
  }
  if (parameter <= knots[degree]) {
    return degree;
  }

  size_t low  = degree;
  size_t high = n;
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
computeBasisFunctions(size_t span, FT parameter, unsigned int degree,
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
    for (unsigned int r = 0; r < j; ++r) {
      FT temp  = basis[r] / (right[r + 1] + left[j - r]);
      basis[r] = saved + right[r + 1] * temp;
      saved    = left[j - r] * temp;
    }
    basis[j] = saved;
  }

  return basis;
}

static auto
computeBasisDerivatives(size_t span, FT parameter, unsigned int degree,
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

    for (unsigned int r = 0; r < j; ++r) {
      ndu[j][r] = right[r + 1] + left[j - r];
      FT temp   = ndu[r][j - 1] / ndu[j][r];
      ndu[r][j] = saved + right[r + 1] * temp;
      saved     = left[j - r] * temp;
    }

    ndu[j][j] = saved;
  }

  for (unsigned int j = 0; j <= degree; ++j) {
    ders[0][j] = ndu[j][degree];
  }

  for (int r = 0; r <= static_cast<int>(degree); ++r) {
    int s1 = 0, s2 = 1;
    a[0][0] = FT(1.0);

    for (unsigned int k = 1; k <= maxDerivative; ++k) {
      FT  d  = FT(0.0);
      int rk = r - static_cast<int>(k);
      int pk = static_cast<int>(degree) - static_cast<int>(k);

      if (r >= static_cast<int>(k)) {
        a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
        d        = a[s2][0] * ndu[rk][pk];
      }

      int j1 = (rk >= -1) ? 1 : -rk;
      int j2 = (r - 1 <= pk) ? static_cast<int>(k) - 1
                             : static_cast<int>(degree) - r;

      for (int j = j1; j <= j2; ++j) {
        a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
        d += a[s2][j] * ndu[rk + j][pk];
      }

      if (r <= pk) {
        a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
        d += a[s2][k] * ndu[r][pk];
      }

      ders[k][r] = d;

      // Switch rows
      std::swap(s1, s2);
    }
  }

  // Multiply through by the correct factors
  int r = static_cast<int>(degree);
  for (unsigned int k = 1; k <= maxDerivative; ++k) {
    for (unsigned int j = 0; j <= degree; ++j) {
      ders[k][j] *= FT(r);
    }
    r *= (static_cast<int>(degree) - static_cast<int>(k));
  }

  return ders;
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

  size_t n = points.size();
  if (n < degree + 1) {
    BOOST_THROW_EXCEPTION(
        Exception("Not enough points for clamped interpolation"));
  }

  // For clamped interpolation, we solve the system N * P = Q
  // where N is the basis function matrix, P are control points, Q are data points
  // The clamped condition ensures the curve passes exactly through endpoints

  matrix<double> N(n, n);

  // Fill basis function matrix
  for (size_t i = 0; i < n; ++i) {
    FT parameter = parameters[i];

    // Find knot span
    size_t span = findKnotSpan(parameter, degree, knots);

    // Compute basis functions at this parameter
    auto basis = computeBasisFunctions(span, parameter, degree, knots);

    // Initialize row to zero
    std::fill(N.data().begin() + i * n, N.data().begin() + (i + 1) * n, 0.0);

    // Fill non-zero entries in this row
    size_t baseIdx = span - degree;
    for (unsigned int j = 0; j <= degree && baseIdx + j < n; ++j) {
      N(i, baseIdx + j) = CGAL::to_double(basis[j]);
    }
  }

  // Solve the system for each coordinate dimension
  std::vector<Point> controlPoints;
  controlPoints.reserve(n);

  // Determine coordinate dimensions
  bool is3D = points.front().is3D();
  bool isMeasured = points.front().isMeasured();

  try {
    // Solve system using LU decomposition for each coordinate
    permutation_matrix<std::size_t> pm(n);
    matrix<double> NCopy;

    // Solve for X coordinates
    vector<double> qx(n);
    for (size_t i = 0; i < n; ++i) {
      qx(i) = CGAL::to_double(points[i].x());
    }
    
    NCopy = N;
    lu_factorize(NCopy, pm);
    lu_substitute(NCopy, pm, qx);
    vector<double> px = qx;

    // Solve for Y coordinates
    vector<double> qy(n);
    for (size_t i = 0; i < n; ++i) {
      qy(i) = CGAL::to_double(points[i].y());
    }
    
    NCopy = N;
    lu_factorize(NCopy, pm);
    lu_substitute(NCopy, pm, qy);
    vector<double> py = qy;

    // Solve for Z coordinates if 3D
    vector<double> pz(n);
    if (is3D) {
      vector<double> qz(n);
      for (size_t i = 0; i < n; ++i) {
        qz(i) = CGAL::to_double(points[i].z());
      }
      
      NCopy = N;
      lu_factorize(NCopy, pm);
      lu_substitute(NCopy, pm, qz);
      pz = qz;
    }

    // Solve for M coordinates if measured
    vector<double> pm_coord(n);
    if (isMeasured) {
      vector<double> qm(n);
      for (size_t i = 0; i < n; ++i) {
        qm(i) = CGAL::to_double(points[i].m());
      }
      
      NCopy = N;
      lu_factorize(NCopy, pm);
      lu_substitute(NCopy, pm, qm);
      pm_coord = qm;
    }

    // Build control points
    for (size_t i = 0; i < n; ++i) {
      FT x = FT(px(i));
      FT y = FT(py(i));
      
      if (is3D && isMeasured) {
        FT z = FT(pz(i));
        double m = pm_coord(i);
        controlPoints.emplace_back(x, y, z, m);  // Uses (FT, FT, FT, double) constructor
      } else if (is3D) {
        FT z = FT(pz(i));
        controlPoints.emplace_back(x, y, z);
      } else if (isMeasured) {
        double x_d = CGAL::to_double(x);
        double y_d = CGAL::to_double(y);
        double m = pm_coord(i);
        controlPoints.emplace_back(x_d, y_d, 0.0, m, COORDINATE_XYM);  // Use explicit constructor with CoordinateType
      } else {
        controlPoints.emplace_back(x, y);
      }
    }

  } catch (const std::exception &e) {
    BOOST_THROW_EXCEPTION(
        Exception("Failed to solve clamped interpolation system: " + std::string(e.what())));
  }

  return controlPoints;
}

auto
NURBSCurve::interpolateNaturalCurve(const std::vector<Point>     &points,
                                    const std::vector<Parameter> &parameters,
                                    unsigned int                  degree,
                                    const std::vector<Knot>      &knots)
    -> std::vector<Point>
{
  using namespace detail::ublas;

  size_t n = points.size();
  if (n < degree + 1) {
    BOOST_THROW_EXCEPTION(
        Exception("Not enough points for natural interpolation"));
  }

  // Set up the interpolation matrix system N * P = Q
  // where N is the basis function matrix, P are control points, Q are data
  // points

  matrix<double> N(n, n);

  // Fill basis function matrix
  for (size_t i = 0; i < n; ++i) {
    FT parameter = parameters[i];

    // Find knot span
    size_t span = findKnotSpan(parameter, degree, knots);

    // Compute basis functions
    auto basis = computeBasisFunctions(span, parameter, degree, knots);

    // Fill row of matrix
    std::fill(N.data().begin() + i * n, N.data().begin() + (i + 1) * n, 0.0);

    size_t baseIdx = span - degree;
    for (unsigned int j = 0; j <= degree && baseIdx + j < n; ++j) {
      N(i, baseIdx + j) = CGAL::to_double(basis[j]);
    }
  }

  // Add natural end conditions (minimize curvature at ends)
  if (degree >= 2) {
    // Replace first and last equations with curvature minimization
    // Second derivative at ends should be minimal

    // For simplicity, use zero second derivative conditions
    // This is a simplified natural condition - full implementation would
    // minimize integral of curvature

    std::fill(N.data().begin(), N.data().begin() + n, 0.0);
    std::fill(N.data().begin() + (n - 1) * n, N.data().begin() + n * n, 0.0);

    // Set up second derivative = 0 at start
    FT     startParam = parameters.front();
    size_t startSpan  = findKnotSpan(startParam, degree, knots);
    auto   startBasis2nd =
        computeBasisDerivatives(startSpan, startParam, degree, knots, 2);

    size_t startBaseIdx = startSpan - degree;
    for (unsigned int j = 0; j <= degree && startBaseIdx + j < n; ++j) {
      N(0, startBaseIdx + j) = CGAL::to_double(startBasis2nd[2][j]);
    }

    // Set up second derivative = 0 at end
    FT     endParam = parameters.back();
    size_t endSpan  = findKnotSpan(endParam, degree, knots);
    auto   endBasis2nd =
        computeBasisDerivatives(endSpan, endParam, degree, knots, 2);

    size_t endBaseIdx = endSpan - degree;
    for (unsigned int j = 0; j <= degree && endBaseIdx + j < n; ++j) {
      N(n - 1, endBaseIdx + j) = CGAL::to_double(endBasis2nd[2][j]);
    }
  }

  // Solve for each coordinate dimension
  std::vector<Point> controlPoints;
  controlPoints.reserve(n);

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
  vector<double> qx(n), px(n);
  for (size_t i = 0; i < n; ++i) {
    qx(i) = CGAL::to_double(points[i].x());
    if (degree >= 2 && (i == 0 || i == n - 1)) {
      qx(i) = 0.0; // Natural end condition: zero second derivative
    }
  }

  // Solve N * px = qx
  permutation_matrix<size_t> pm(n);
  matrix<double>             NCopy = N;
  lu_factorize(NCopy, pm);
  lu_substitute(NCopy, pm, qx);
  px = qx; // qx now contains solution

  // Solve for Y coordinates
  vector<double> qy(n), py(n);
  for (size_t i = 0; i < n; ++i) {
    qy(i) = CGAL::to_double(points[i].y());
    if (degree >= 2 && (i == 0 || i == n - 1)) {
      qy(i) = 0.0; // Natural end condition
    }
  }

  NCopy = N;
  lu_factorize(NCopy, pm);
  lu_substitute(NCopy, pm, qy);
  py = qy;

  // Solve for Z coordinates if 3D
  vector<double> pz(n);
  if (points[0].is3D()) {
    vector<double> qz(n);
    for (size_t i = 0; i < n; ++i) {
      qz(i) = CGAL::to_double(points[i].z());
      if (degree >= 2 && (i == 0 || i == n - 1)) {
        qz(i) = 0.0;
      }
    }

    NCopy = N;
    lu_factorize(NCopy, pm);
    lu_substitute(NCopy, pm, qz);
    pz = qz;
  }

  // Solve for M coordinates if measured
  vector<double> pm_coord(n);
  if (points[0].isMeasured()) {
    vector<double> qm(n);
    for (size_t i = 0; i < n; ++i) {
      qm(i) = CGAL::to_double(points[i].m());
      if (degree >= 2 && (i == 0 || i == n - 1)) {
        qm(i) = 0.0;
      }
    }

    NCopy = N;
    lu_factorize(NCopy, pm);
    lu_substitute(NCopy, pm, qm);
    pm_coord = qm;
  }

  // Build control points
  for (size_t i = 0; i < n; ++i) {
    FT     x = FT(px(i));
    FT     y = FT(py(i));
    FT     z = points[0].is3D() ? FT(pz(i)) : FT(0);
    double m = points[0].isMeasured() ? pm_coord(i) : NaN();

    controlPoints.emplace_back(x, y, z, m, coordType);
  }

  return controlPoints;
}

auto
NURBSCurve::interpolatePeriodicCurve(const std::vector<Point> &points,
                                     unsigned int degree, KnotMethod knotMethod)
    -> std::unique_ptr<NURBSCurve>
{
  if (points.size() < degree + 1) {
    BOOST_THROW_EXCEPTION(
        Exception("Not enough points for periodic interpolation"));
  }

  // Check if curve is actually closed
  FT tolerance = FT(1e-10);
  if (algorithm::distance(points.front(), points.back()) > tolerance) {
    BOOST_THROW_EXCEPTION(
        Exception("Points must form closed curve for periodic interpolation"));
  }

  // For periodic curves, we need to wrap points and create uniform knot spacing
  std::vector<Point> periodicPoints = points;

  // Remove duplicate end point for internal processing
  if (periodicPoints.size() > 1) {
    periodicPoints.pop_back();
  }

  size_t n = periodicPoints.size();

  // Create periodic parameter vector
  std::vector<Parameter> parameters;
  parameters.reserve(n);

  for (size_t i = 0; i < n; ++i) {
    parameters.push_back(FT(i) / FT(n));
  }

  // Create periodic knot vector
  std::vector<Knot> knots;
  size_t            numKnots = n + degree + 1;
  knots.reserve(numKnots);

  // Periodic knots wrap around
  for (size_t i = 0; i < numKnots; ++i) {
    int knotIndex = static_cast<int>(i) - static_cast<int>(degree);
    knots.push_back(FT(knotIndex) / FT(n));
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
NURBSCurve::generateApproximationKnotVector(const std::vector<Parameter> &parameters,
                                           unsigned int degree, size_t numControlPoints)
    -> std::vector<Knot>
{
  // Generate knot vector for approximation (Piegl & Tiller, "The NURBS Book", Algorithm A9.1)
  size_t n = parameters.size() - 1;  // Number of data points - 1
  size_t m = numControlPoints - 1;   // Number of control points - 1

  std::vector<Knot> knots;
  knots.reserve(m + degree + 2);     // Correct size: m + p + 1 + 1

  // First degree+1 knots are 0
  for (unsigned int i = 0; i <= degree; ++i) {
    knots.push_back(FT(0));
  }

  // Internal knots (Piegl & Tiller averaging method)
  // We need exactly (m - degree + 1) internal knots for proper size
  for (size_t j = 1; j <= m - degree + 1; ++j) {
    FT sum = FT(0);
    for (unsigned int k = 0; k < degree; ++k) {
      // Use direct indexing into parameter array (Piegl & Tiller, Equation 9.8)
      size_t dataIdx = j + k;
      if (dataIdx < parameters.size()) {
        sum += parameters[dataIdx];
      }
    }
    knots.push_back(sum / FT(degree));
  }

  // Last degree+1 knots are 1
  for (unsigned int i = 0; i <= degree; ++i) {
    knots.push_back(FT(1));
  }

  return knots;
}

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
  numControlPoints = std::max(numControlPoints, minControlPoints);
  
  if (numControlPoints >= points.size()) {
    // If we need as many control points as data points, use interpolation instead
    return interpolateCurve(points, degree, KnotMethod::CHORD_LENGTH, EndCondition::CLAMPED);
  }

  // NURBS Curve Approximation Algorithm (Piegl & Tiller, "The NURBS Book", Chapter 9)
  // This implements true least-squares approximation, NOT sampling

  using namespace detail::ublas;
  
  size_t n = points.size() - 1;  // Number of data points - 1 (index range)
  size_t m = numControlPoints - 1;  // Number of control points - 1

  // Step 1: Compute parameter values for data points (chord length parameterization)
  std::vector<Parameter> parameters = computeParameters(points, KnotMethod::CHORD_LENGTH);

  // Step 2: Generate knot vector for approximation
  std::vector<Knot> knots = generateApproximationKnotVector(parameters, degree, m);

  // Step 3: Set up the least-squares system N^T * N * P = N^T * Q
  // where N is the basis function matrix, P are control points, Q are data points

  matrix<double> N(n + 1, m + 1);  // Basis function matrix
  
  // Fill basis function matrix
  for (size_t i = 0; i <= n; ++i) {
    Parameter t = parameters[i];
    
    // Find knot span
    size_t span = findKnotSpan(t, degree, knots);
    
    // Compute non-zero basis functions
    auto basis = computeBasisFunctions(span, t, degree, knots);
    
    // Initialize row to zero
    for (size_t j = 0; j <= m; ++j) {
      N(i, j) = 0.0;
    }
    
    // Fill non-zero basis functions
    size_t baseIdx = span - degree;
    for (unsigned int k = 0; k <= degree && baseIdx + k <= m; ++k) {
      N(i, baseIdx + k) = CGAL::to_double(basis[k]);
    }
  }

  // Step 4: Compute N^T * N (normal matrix)
  matrix<double> NTN = prod(trans(N), N);

  // Step 5: Solve for control points in each coordinate dimension
  std::vector<Point> controlPoints;
  controlPoints.reserve(m + 1);

  // Determine coordinate dimensions
  bool is3D = points.front().is3D();
  bool isMeasured = points.front().isMeasured();

  try {
    // Solve for X coordinates
    vector<double> qx(n + 1);
    for (size_t i = 0; i <= n; ++i) {
      qx(i) = CGAL::to_double(points[i].x());
    }
    vector<double> ntqx = prod(trans(N), qx);  // N^T * Q_x
    
    // Solve NTN * px = ntqx
    permutation_matrix<std::size_t> pm(m + 1);
    matrix<double> NTN_copy = NTN;
    lu_factorize(NTN_copy, pm);
    lu_substitute(NTN_copy, pm, ntqx);
    vector<double> px = ntqx;

    // Solve for Y coordinates
    vector<double> qy(n + 1);
    for (size_t i = 0; i <= n; ++i) {
      qy(i) = CGAL::to_double(points[i].y());
    }
    vector<double> ntqy = prod(trans(N), qy);  // N^T * Q_y
    
    NTN_copy = NTN;
    lu_factorize(NTN_copy, pm);
    lu_substitute(NTN_copy, pm, ntqy);
    vector<double> py = ntqy;

    // Solve for Z coordinates if 3D
    vector<double> pz(m + 1, 0.0);
    if (is3D) {
      vector<double> qz(n + 1);
      for (size_t i = 0; i <= n; ++i) {
        qz(i) = CGAL::to_double(points[i].z());
      }
      vector<double> ntqz = prod(trans(N), qz);  // N^T * Q_z
      
      NTN_copy = NTN;
      lu_factorize(NTN_copy, pm);
      lu_substitute(NTN_copy, pm, ntqz);
      pz = ntqz;
    }

    // Solve for M coordinates if measured
    vector<double> pm_coord(m + 1, 0.0);
    if (isMeasured) {
      vector<double> qm(n + 1);
      for (size_t i = 0; i <= n; ++i) {
        qm(i) = CGAL::to_double(points[i].m());
      }
      vector<double> ntqm = prod(trans(N), qm);  // N^T * Q_m
      
      NTN_copy = NTN;
      lu_factorize(NTN_copy, pm);
      lu_substitute(NTN_copy, pm, ntqm);
      pm_coord = ntqm;
    }

    // Build control points from solved coordinates
    for (size_t i = 0; i <= m; ++i) {
      FT x = FT(px(i));
      FT y = FT(py(i));
      
      if (is3D && isMeasured) {
        FT z = FT(pz(i));
        double m_val = pm_coord(i);
        controlPoints.emplace_back(x, y, z, m_val);
      } else if (is3D) {
        FT z = FT(pz(i));
        controlPoints.emplace_back(x, y, z);
      } else if (isMeasured) {
        double x_d = CGAL::to_double(x);
        double y_d = CGAL::to_double(y);
        double m_val = pm_coord(i);
        controlPoints.emplace_back(x_d, y_d, 0.0, m_val, COORDINATE_XYM);
      } else {
        controlPoints.emplace_back(x, y);
      }
    }

  } catch (const std::exception &e) {
    BOOST_THROW_EXCEPTION(
        Exception("Failed to solve least-squares approximation system: " + std::string(e.what())));
  }

  // Create the approximating curve
  auto curve = std::make_unique<NURBSCurve>(
      controlPoints, std::vector<FT>(controlPoints.size(), FT(1)), 
      degree, knots);
  
  // Store fit information
  curve->_fitPoints = points;
  curve->_fitTolerance = tolerance;

  return curve;
}

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
    return approximateCurve(points, degree, tolerance, maxControlPoints);

  default:
    BOOST_THROW_EXCEPTION(Exception("Unknown FitMethod"));
  }
}

//-- Geometry interface implementation

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

//-- Core curve evaluation and derivatives

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
    return Point(FT(0), FT(0), FT(0), 0.0);
  } else if (is3D()) {
    return Point(FT(0), FT(0), FT(0));
  } else if (isMeasured()) {
    Point point(FT(0), FT(0));
    point.setM(0.0);
    return point;
  } else {
    return Point(FT(0), FT(0));
  }
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

  if (magnitude < FT(1e-12)) {
    BOOST_THROW_EXCEPTION(
        Exception("Cannot compute tangent at singular point"));
  }

  return Point(firstDerivative.x() / magnitude, firstDerivative.y() / magnitude,
               firstDerivative.is3D() ? firstDerivative.z() / magnitude
                                      : FT(0));
}

auto
NURBSCurve::normal(Parameter parameter) const -> Point
{
  if (is3D()) {
    BOOST_THROW_EXCEPTION(
        Exception("Normal vector not uniquely defined for 3D curves"));
  }

  Point tangentVec = tangent(parameter);
  return Point(-tangentVec.y(), tangentVec.x());
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

  if (magnitude < FT(1e-12)) {
    BOOST_THROW_EXCEPTION(
        Exception("Cannot compute binormal at inflection point"));
  }

  return Point(crossX / magnitude, crossY / magnitude, crossZ / magnitude);
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

    if (firstMagnitude < FT(1e-12)) {
      return FT(0);
    }

    return crossMagnitude / (firstMagnitude * firstMagnitude * firstMagnitude);
  } else {
    FT numerator =
        firstDeriv.x() * secondDeriv.y() - firstDeriv.y() * secondDeriv.x();
    FT denominator =
        firstDeriv.x() * firstDeriv.x() + firstDeriv.y() * firstDeriv.y();

    if (denominator < FT(1e-12)) {
      return FT(0);
    }

    return numerator / std::pow(CGAL::to_double(denominator), 1.5);
  }
}

auto
NURBSCurve::torsion(Parameter parameter) const -> FT
{
  if (!is3D()) {
    return FT(0);
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

  if (crossMagnitudeSquared < FT(1e-12)) {
    return FT(0);
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
  } else {
    Point normalVec   = normal(parameter);
    Point binormalVec = Point(0, 0, 1);
    return std::make_tuple(tangentVec, normalVec, binormalVec);
  }
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
NURBSCurve::toLineStringArcLength(unsigned int numSegments) const
    -> std::unique_ptr<LineString>
{
  auto lineString = std::make_unique<LineString>();
  if (_controlPoints.empty()) {
    return lineString;
  }

  if (numSegments == 0) {
    numSegments = 32;
  }

  // For performance reasons, limit the number of segments for arc-length sampling
  if (numSegments > 100) {
    // Fall back to parameter-based sampling for very high segment counts
    return toLineString(numSegments);
  }

  // Calculate total arc length
  FT totalLength = length();
  if (totalLength <= FT(1e-12)) {
    // Degenerate curve, fall back to parameter sampling
    return toLineString(numSegments);
  }

  auto bounds = parameterBounds();
  
  // Build a lookup table of parameter values vs cumulative arc lengths
  // Use more samples than requested to build accurate inverse mapping
  const unsigned int lookupSamples = std::max(numSegments * 2, 64u);
  std::vector<Parameter> lookupParams;
  std::vector<FT> lookupLengths;
  lookupParams.reserve(lookupSamples + 1);
  lookupLengths.reserve(lookupSamples + 1);
  
  FT paramStep = (bounds.second - bounds.first) / FT(lookupSamples);
  lookupParams.push_back(bounds.first);
  lookupLengths.push_back(FT(0));
  
  for (unsigned int i = 1; i <= lookupSamples; ++i) {
    Parameter param = bounds.first + FT(i) * paramStep;
    FT cumulativeLength = computeArcLength(bounds.first, param, FT(1e-6));
    lookupParams.push_back(param);
    lookupLengths.push_back(cumulativeLength);
  }

  // Now sample points uniformly in arc-length space using interpolation
  lineString->addPoint(evaluate(bounds.first));
  
  for (unsigned int segIdx = 1; segIdx <= numSegments; ++segIdx) {
    FT targetArcLength = (FT(segIdx) / FT(numSegments)) * totalLength;
    
    // Find the parameter by linear interpolation in the lookup table
    Parameter param = bounds.first;
    
    // Find bracketing indices in lookup table
    for (size_t j = 0; j < lookupLengths.size() - 1; ++j) {
      if (targetArcLength >= lookupLengths[j] && targetArcLength <= lookupLengths[j + 1]) {
        // Linear interpolation
        FT lengthRange = lookupLengths[j + 1] - lookupLengths[j];
        if (lengthRange > FT(1e-12)) {
          FT ratio = (targetArcLength - lookupLengths[j]) / lengthRange;
          param = lookupParams[j] + ratio * (lookupParams[j + 1] - lookupParams[j]);
        } else {
          param = lookupParams[j];
        }
        break;
      }
    }
    
    lineString->addPoint(evaluate(param));
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
        parameters.insert(parameters.begin() + paramIdx + 1, midParam);
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
        FT(1e-10)) {
      return false;
    }
  }

  return true;
}

auto
NURBSCurve::isPlanar(std::vector<FT> *plane) const -> bool
{
  if (!is3D() || _controlPoints.size() < 4) {
    if (plane) {
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

        if (normalMagnitude > FT(1e-10)) {
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
            if (distance > FT(1e-10)) {
              allPointsOnPlane = false;
              break;
            }
          }

          if (allPointsOnPlane) {
            if (plane) {
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

  if (dirMagnitude < FT(1e-10)) {
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

    if (vecMagnitude > FT(1e-10)) {
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

      if (crossMagnitude > FT(1e-10)) {
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

  if (from < FT(0))
    from = bounds.first;
  if (to < FT(0))
    to = bounds.second;

  if (from >= to) {
    return FT(0);
  }

  tolerance = CGAL::abs(tolerance);
  if (tolerance <= FT(0)) {
    tolerance = FT(1e-6);
  }

  return computeArcLength(from, to, tolerance);
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

  if (totalLength <= FT(1e-12)) {
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
    Parameter param        = parameterAtLength(targetLength, FT(1e-8));
    newControlPoints.push_back(evaluate(param));
  }

  // Use a reasonable degree for the reparameterized curve
  unsigned int newDegree =
      std::min(3u, static_cast<unsigned int>(newControlPoints.size() - 1));

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
  from        = std::max(from, bounds.first);
  to          = std::min(to, bounds.second);

  if (from >= to) {
    return std::make_unique<NURBSCurve>();
  }

  unsigned int numSamples =
      std::max(16u, static_cast<unsigned int>(_controlPoints.size()));

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
  } else {
    return std::make_unique<NURBSCurve>(subPoints, subWeights, subDegree,
                                        subKnots);
  }
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
    FT maxKnot = _knotVector.back();

    for (auto knotIter = _knotVector.rbegin(); knotIter != _knotVector.rend();
         ++knotIter) {
      reversedKnots.push_back(maxKnot - *knotIter);
    }
  }

  return std::make_unique<NURBSCurve>(reversedPoints, reversedWeights, _degree,
                                      reversedKnots);
}

auto
NURBSCurve::join(const Curve &other, Continuity continuity, FT tolerance) const
    -> std::unique_ptr<Curve>
{
  const auto *otherNurbs = dynamic_cast<const NURBSCurve *>(&other);
  if (!otherNurbs) {
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

  std::vector<Point> joinedPoints = _controlPoints;
  const auto        &otherPoints  = otherNurbs->controlPoints();

  joinedPoints.insert(joinedPoints.end(), otherPoints.begin() + 1,
                      otherPoints.end());

  unsigned int joinedDegree = std::max(_degree, otherNurbs->degree());
  joinedDegree              = std::min(joinedDegree,
                                       static_cast<unsigned int>(joinedPoints.size() - 1));

  return std::make_unique<NURBSCurve>(joinedPoints, joinedDegree);
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

  if (outParameter) {
    *outParameter = parameter;
  }

  return closestPt;
}

//-- NURBS-specific data access and manipulation

void
NURBSCurve::setControlPoints(const std::vector<Point> &controlPoints)
{
  _controlPoints = controlPoints;

  if (!_weights.empty() && _weights.size() != controlPoints.size()) {
    _weights.resize(controlPoints.size(), FT(1));
  }

  if (!_knotVector.empty() && !controlPoints.empty()) {
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
    return FT(1);
  }
  if (index >= _weights.size()) {
    BOOST_THROW_EXCEPTION(Exception("Weight index out of bounds"));
  }
  return _weights[index];
}

void
NURBSCurve::setWeight(size_t index, FT weight)
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
  const FT tolerance   = FT(1e-10);

  for (const auto &weight : _weights) {
    if (CGAL::abs(weight - firstWeight) > tolerance) {
      return true;
    }
  }

  return false;
}

auto
NURBSCurve::isBezier() const -> bool
{
  if (_knotVector.empty() || _controlPoints.empty()) {
    return false;
  }

  size_t expectedSize = 2 * (_degree + 1);
  if (_knotVector.size() != expectedSize) {
    return false;
  }

  const FT tolerance = FT(1e-10);

  for (unsigned int idx = 0; idx <= _degree; ++idx) {
    if (CGAL::abs(_knotVector[idx] - _knotVector[0]) > tolerance) {
      return false;
    }
  }

  for (unsigned int idx = _degree + 1; idx < _knotVector.size(); ++idx) {
    if (CGAL::abs(_knotVector[idx] - _knotVector.back()) > tolerance) {
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
NURBSCurve::knotMultiplicity(Knot value, FT tolerance) const -> unsigned int
{
  unsigned int multiplicity = 0;
  for (const auto &knot : _knotVector) {
    if (CGAL::abs(knot - value) <= tolerance) {
      ++multiplicity;
    }
  }
  return multiplicity;
}

//-- Advanced NURBS operations

auto
NURBSCurve::insertKnot(Knot parameter, unsigned int times) const
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
NURBSCurve::reduceDegree(FT tolerance) const -> std::unique_ptr<NURBSCurve>
{
  BOOST_THROW_EXCEPTION(
      Exception("reduceDegree: NURBS degree reduction not yet implemented"));
}

auto
NURBSCurve::removeKnots(FT tolerance) const -> std::unique_ptr<NURBSCurve>
{
  BOOST_THROW_EXCEPTION(
      Exception("removeKnots: NURBS knot removal not yet implemented"));
}

//-- Validation and utility methods

auto
NURBSCurve::validateData() const -> std::pair<bool, std::string>
{
  return validateNURBSData(_controlPoints, _weights, _degree, _knotVector);
}

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

  if (!knots.empty() && controlPoints.size() > 0) {
    FT startParam = knots[degree];
    FT endParam   = knots[controlPoints.size()];

    if (startParam > endParam) {
      return {false, "Invalid parameter bounds"};
    }
  }

  return {true, ""};
}

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

//-- Protected helper methods

auto
NURBSCurve::findSpan(Parameter parameter) const -> size_t
{
  // Binary search algorithm to find knot span containing parameter
  // Based on Algorithm A2.1 from "The NURBS Book" by Piegl & Tiller
  // Returns index i such that t_i <= parameter < t_{i+1}

  size_t numControlPoints = _controlPoints.size();

  // Handle boundary cases - parameter at or beyond curve end
  if (parameter >= _knotVector[numControlPoints]) {
    return numControlPoints - 1;
  }

  // Handle boundary cases - parameter at or before curve start
  if (parameter <= _knotVector[_degree]) {
    return _degree;
  }

  // Binary search for the knot span
  size_t low  = _degree;
  size_t high = numControlPoints;
  size_t mid  = (low + high) / 2;

  while (parameter < _knotVector[mid] || parameter >= _knotVector[mid + 1]) {
    if (parameter < _knotVector[mid]) {
      high = mid;
    } else {
      low = mid;
    }
    mid = (low + high) / 2;
  }

  return mid;
}

auto
NURBSCurve::deBoorRational(size_t span, Parameter parameter) const -> Point
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

    HomogeneousPoint(const Point &point, FT pointWeight)
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

      if (CGAL::abs(denominator) > FT(1e-15)) {
        alpha = (parameter - _knotVector[knotIdx]) / denominator;
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

  if (CGAL::abs(result.weight) < FT(1e-15)) {
    BOOST_THROW_EXCEPTION(Exception("Division by zero in NURBS evaluation"));
  }

  FT invWeight = FT(1) / result.weight;

  if (is3D() && isMeasured()) {
    return Point(result.weightedX * invWeight, result.weightedY * invWeight,
                 result.weightedZ * invWeight, result.measure);
  } else if (is3D()) {
    return Point(result.weightedX * invWeight, result.weightedY * invWeight,
                 result.weightedZ * invWeight);
  } else if (isMeasured()) {
    Point point(result.weightedX * invWeight, result.weightedY * invWeight);
    point.setM(result.measure);
    return point;
  } else {
    return Point(result.weightedX * invWeight, result.weightedY * invWeight);
  }
}

auto
NURBSCurve::computeDerivatives(Parameter parameter, unsigned int maxOrder) const
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
      const FT epsilon = FT(1e-8);
      auto     bounds  = parameterBounds();

      for (unsigned int orderIdx = 1; orderIdx <= maxOrder; ++orderIdx) {
        FT param1     = std::max(bounds.first, parameter - epsilon);
        FT param2     = std::min(bounds.second, parameter + epsilon);
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
      knots.push_back(FT(pointIdx));
    }
    knots.push_back(FT(numPoints));
    return knots;
  }

  knots.reserve(numKnots);

  if (method == KnotMethod::UNIFORM) {
    for (unsigned int degreeIdx = 0; degreeIdx <= degree; ++degreeIdx) {
      knots.push_back(FT(0));
    }

    size_t numInternal = numKnots - 2 * (degree + 1);
    for (size_t internalIdx = 1; internalIdx <= numInternal; ++internalIdx) {
      knots.push_back(FT(internalIdx) / FT(numInternal + 1));
    }

    for (unsigned int degreeIdx = 0; degreeIdx <= degree; ++degreeIdx) {
      knots.push_back(FT(1));
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
    parameters.push_back(FT(0));
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

  for (const auto &point : _controlPoints) {
    if (point.is3D() != first3D || point.isMeasured() != firstMeasured) {
      return false;
    }
  }

  return true;
}

auto
NURBSCurve::computeArcLength(Parameter startParam, Parameter endParam,
                             FT tolerance) const -> FT
{
  if (startParam >= endParam) {
    return FT(0);
  }

  // Use adaptive Simpson's rule for arc length integration
  // This integrates ||C'(t)|| dt over [startParam, endParam]
  return adaptiveSimpsonArcLength(startParam, endParam, tolerance);
}

//-- Private helper methods for adaptive Simpson integration

auto
NURBSCurve::speedFunction(Parameter t) const -> FT
{
  // Compute ||C'(t)|| - the magnitude of the derivative
  Point derivative = this->derivative(t, 1);
  
  if (is3D()) {
    FT dx = derivative.x();
    FT dy = derivative.y(); 
    FT dz = derivative.z();
    // Use CGAL::to_double for sqrt calculation, then convert back to FT
    double magnitude = std::sqrt(CGAL::to_double(dx*dx + dy*dy + dz*dz));
    return FT(magnitude);
  } else {
    FT dx = derivative.x();
    FT dy = derivative.y();
    // Use CGAL::to_double for sqrt calculation, then convert back to FT
    double magnitude = std::sqrt(CGAL::to_double(dx*dx + dy*dy));
    return FT(magnitude);
  }
}

auto
NURBSCurve::simpsonRule(Parameter a, Parameter b) const -> FT
{
  // Basic Simpson's rule: (b-a)/6 * [f(a) + 4*f((a+b)/2) + f(b)]
  Parameter mid = (a + b) / FT(2);
  FT fa = speedFunction(a);
  FT fmid = speedFunction(mid);
  FT fb = speedFunction(b);
  
  return (b - a) / FT(6) * (fa + FT(4) * fmid + fb);
}

auto
NURBSCurve::adaptiveSimpsonArcLength(Parameter a, Parameter b, FT tolerance,
                                     unsigned int maxDepth) const -> FT
{
  if (maxDepth == 0) {
    // Fallback to basic Simpson's rule when max recursion reached
    return simpsonRule(a, b);
  }

  Parameter mid = (a + b) / FT(2);
  
  // Compute Simpson's rule for whole interval and two halves
  FT wholeInterval = simpsonRule(a, b);
  FT leftHalf = simpsonRule(a, mid);
  FT rightHalf = simpsonRule(mid, b);
  FT twoHalves = leftHalf + rightHalf;
  
  // Error estimate: |S(a,b) - S(a,m) - S(m,b)| / 15
  FT error = CGAL::abs(wholeInterval - twoHalves) / FT(15);
  
  if (error <= tolerance) {
    // Richardson extrapolation: more accurate result
    return twoHalves + (twoHalves - wholeInterval) / FT(15);
  } else {
    // Recursively subdivide with tighter tolerance
    FT halfTolerance = tolerance / FT(2);
    return adaptiveSimpsonArcLength(a, mid, halfTolerance, maxDepth - 1) +
           adaptiveSimpsonArcLength(mid, b, halfTolerance, maxDepth - 1);
  }
}

auto
NURBSCurve::findParameterByArcLength(FT targetLength, FT tolerance) const
    -> Parameter
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
    if (high - low < FT(1e-12)) {
      break;
    }
  }

  return (low + high) / FT(2);
}

auto
NURBSCurve::projectPointToCurve(const Point &point, FT tolerance,
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

    Point diff = Point(point.x() - curvePoint.x(), point.y() - curvePoint.y(),
                       point.is3D() ? point.z() - curvePoint.z() : FT(0));

    FT numerator = diff.x() * firstDeriv.x() + diff.y() * firstDeriv.y() +
                   (diff.is3D() ? diff.z() * firstDeriv.z() : FT(0));

    FT denominator =
        firstDeriv.x() * firstDeriv.x() + firstDeriv.y() * firstDeriv.y() +
        (firstDeriv.is3D() ? firstDeriv.z() * firstDeriv.z() : FT(0));

    if (denominator < FT(1e-12)) {
      break;
    }

    FT paramUpdate = numerator / denominator;

    if (CGAL::abs(paramUpdate) <= tolerance) {
      break;
    }

    bestParam += paramUpdate;
    bestParam = std::max(bounds.first, std::min(bounds.second, bestParam));
  }

  return std::make_pair(evaluate(bestParam), bestParam);
}

void
NURBSCurve::reserveCapacity(size_t expectedSize)
{
  _controlPoints.reserve(expectedSize);
  _weights.reserve(expectedSize);
  _knotVector.reserve(expectedSize + _degree + 1);
  _fitPoints.reserve(expectedSize);
}

void
NURBSCurve::clearOptionalData()
{
  _fitPoints.clear();
  _fitPoints.shrink_to_fit();
  _startTangent.reset();
  _endTangent.reset();
}

} // namespace SFCGAL
