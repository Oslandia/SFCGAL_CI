// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/primitive3d/Primitive.h"
#include "SFCGAL/Exception.h"

namespace SFCGAL {

auto
Primitive::operator=(const Primitive &other) -> Primitive &
{
  if (this != &other) {
    m_parameters         = other.m_parameters;
    m_polyhedral_surface = other.m_polyhedral_surface;
  }
  return *this;
}

void
Primitive::setParameter(const std::string        &name,
                        const PrimitiveParameter &parameter)
{
  auto parameterIt = m_parameters.find(name);

  // check that parameter name exists
  if (parameterIt == m_parameters.end()) {
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("%s does not have a parameter named %s") %
                   primitiveType() % name)
                      .str()));
  }

  // check the type of the parameter
  const bool typeMatches = std::visit(
      [&](auto &&currentParameter) {
        using expectedType = std::decay_t<decltype(currentParameter)>;
        return std::holds_alternative<expectedType>(parameter);
      },
      parameterIt->second);

  if (!typeMatches) {
    BOOST_THROW_EXCEPTION(Exception(
        (boost::format("Wrong type for parameter %s.") % name).str()));
    return;
  }

  // Validate the new value
  validateAndSetParameter(name, parameter);
}

void
Primitive::validateAndSetParameter(const std::string        &name,
                                   const PrimitiveParameter &parameter)
{
  // check if valid
  auto tempParameters     = m_parameters;
  tempParameters.at(name) = parameter;
  validateParameters(tempParameters);

  // assign parameter
  m_parameters.at(name) = parameter;

  // call propagation callback
  onValidatedAndSetParameter(name, parameter);

  invalidateCache();
}

void
Primitive::onValidatedAndSetParameter(const std::string        &name,
                                      const PrimitiveParameter &parameter)
{
  (void)name;
  (void)parameter;
  // nope: does nothing by default
}

auto
Primitive::parameter(const std::string &name) const -> PrimitiveParameter
{
  auto parameterIt = m_parameters.find(name);
  // check that parameter name exists
  if (parameterIt == m_parameters.end()) {
    BOOST_THROW_EXCEPTION(
        Exception((boost::format("%s does not have a parameter named %s") %
                   primitiveType() % name)
                      .str()));
  }

  return parameterIt->second;
}

auto
Primitive::parameters() const
    -> std::unordered_map<std::string, PrimitiveParameter>
{
  return m_parameters;
}

void
Primitive::invalidateCache()
{
  m_polyhedral_surface.reset();
}

auto
operator==(const Primitive &pa, const Primitive &pb) -> bool
{
  return pa.primitiveTypeId() == pb.primitiveTypeId() &&
         pa.parameters() == pb.parameters();
}

} // namespace SFCGAL
