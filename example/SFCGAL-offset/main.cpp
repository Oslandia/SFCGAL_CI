// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/MultiPolygon.h"

#include "SFCGAL/algorithm/offset.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;

int
main()
{
  std::unique_ptr<Geometry> g(
      io::readWkt("MULTIPOINT (0 0,5 6,3 2,7 1,4 1,3 5,2 9)"));

  for (size_t i = 1; i <= 50; i++) {
    std::unique_ptr<Geometry> buffer(algorithm::offset(*g, 0.2 * i));
    std::cout << buffer->asText(5) << std::endl;
  }

  return 0;
}
