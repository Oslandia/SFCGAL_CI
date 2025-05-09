#ifndef SFCGALOP_OPERATIONS_HPP
#define SFCGALOP_OPERATIONS_HPP

#include "../io.hpp"
#include <SFCGAL/Geometry.h>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace Operations {

// Type alias for operation metadata: (name, category, description, help,
// params)
using OperationInfo =
    std::tuple<std::string, std::string, std::string, std::string, std::string>;

auto
execute_operation(const std::string &op_name, const std::string &op_arg,
                  const SFCGAL::Geometry *geom_a,
                  const SFCGAL::Geometry *geom_b)
    -> std::optional<OperationResult>;

void
print_available_operations();
auto
print_operation_help(const char *name) -> bool;
auto
get_all_operations_info() -> std::vector<OperationInfo>;

auto
operation_requires_second_geometry(const std::string &operation_name) -> bool;

} // namespace Operations

// Provide non-namespaced aliases for backward compatibility and convenience.
// These allow calling execute_operation() without the Operations:: prefix.
// If namespace pollution is a concern, use the fully qualified names instead.
using Operations::execute_operation;
using Operations::print_available_operations;
using Operations::print_operation_help;

#endif // SFCGALOP_OPERATIONS_HPP