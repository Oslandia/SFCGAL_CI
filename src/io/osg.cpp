// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/config.h"
#ifdef SFCGAL_WITH_OSG

  #include "SFCGAL/detail/io/OsgFactory.h"
  #include "SFCGAL/io/osg.h"

  #include <osg/Geode>
  #include <osg/Geometry>
  #include <osg/Group>

  #include <osgDB/WriteFile>

namespace SFCGAL {
namespace io {

void
osgWriteFile(const Geometry &g, const std::string &filepath)
{
  SFCGAL::detail::io::OsgFactory factory;
  osg::ref_ptr<osg::Geometry>    osgGeometry = factory.createGeometry(g);
  osg::ref_ptr<osg::Geode>       geode       = new osg::Geode;
  geode->setName(g.geometryType());
  geode->addDrawable(osgGeometry);
  osgDB::writeNodeFile(*geode, filepath);
}

osg::Geometry *
toOsgGeometry(const Geometry &g)
{
  SFCGAL::detail::io::OsgFactory factory;
  return factory.createGeometry(g);
}

} // namespace io
} // namespace SFCGAL

#endif // SFCGAL_WITH_OSG
