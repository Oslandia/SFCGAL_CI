#include "SFCGAL/Curve.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/config.h"

namespace SFCGAL {

auto
Curve::startPoint() const -> Point
{
  auto bounds = parameterBounds();
  return evaluate(bounds.first);
}

auto
Curve::endPoint() const -> Point
{
  auto bounds = parameterBounds();
  return evaluate(bounds.second);
}

} // namespace SFCGAL
