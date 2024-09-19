#ifndef _SFCGAL_ALGORITHM_LEXICOGRAPHICORDER_H_
#define _SFCGAL_ALGORITHM_LEXICOGRAPHICORDER_H_

#include <SFCGAL/Geometry.h>
#include <memory>

namespace SFCGAL {
namespace transform {

std::unique_ptr<Geometry>lexicographicOrder(const Geometry &geometry);

}
} // namespace SFCGAL

#endif
