// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/connection.h"

#include "SFCGAL/Coordinate.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

#include <limits>

namespace SFCGAL::algorithm {

const size_t SurfaceGraph::INVALID_INDEX = std::numeric_limits<size_t>::max();

void
SurfaceGraph::addRing(const LineString &ring, FaceIndex faceIndex)
{
  const size_t numSegments = ring.numSegments();

  for (size_t s = 0; s != numSegments; ++s) { // for each segment
    const Coordinate startCoord = ring.pointN(s).coordinate();
    const Coordinate endCoord =
        ring.pointN((s + 1) % numSegments).coordinate(); // possible
                                                         // optimization: store
                                                         // the index of ring
                                                         // start point instead
                                                         // of finding it
    const auto startFound = _coordinateMap.find(startCoord);
    const auto endFound   = _coordinateMap.find(endCoord);
    BOOST_ASSERT(s + 1 != numSegments ||
                 endFound != _coordinateMap.end()); // ring not closed

    if (startFound != _coordinateMap.end() &&
        endFound != _coordinateMap.end()) {
      // found both end, we look for the edge
      const VertexIndex                         startIndex = startFound->second;
      const VertexIndex                         endIndex   = endFound->second;
      const std::pair<VertexIndex, VertexIndex> edge(startIndex, endIndex);
      const auto foundEdgeWithBadOrientation = _edgeMap.find(edge);

      if (foundEdgeWithBadOrientation != _edgeMap.end()) {
        _isValid = Validity::invalid(
            (boost::format("inconsistent orientation of PolyhedralSurface "
                           "detected at edge %d (%d-%d) of polygon %d") %
             s % edge.first % edge.second % faceIndex)
                .str());
      }

      const std::pair<VertexIndex, VertexIndex> reversedEdge(endIndex,
                                                             startIndex);

      const auto foundEdge = _edgeMap.find(reversedEdge);

      if (foundEdge != _edgeMap.end()) {
        // edit edge
        foundEdge->second.second = faceIndex;
        // we have two faces connected, this is an edge of the graph
        boost::add_edge(foundEdge->second.first, foundEdge->second.second,
                        _graph);
        // std::cerr << "face " << foundEdge->second.first << "->" <<
        // foundEdge->second.second << "\n";
      } else {
        // create edge
        _edgeMap.insert(
            std::make_pair(edge, std::make_pair(faceIndex, INVALID_INDEX)));
        // std::cerr << "face " << faceIndex << " edge " << edge.first << "-" <<
        // edge.second << "\n";
      }
    } else {
      // one end at least is missing, create the edge
      VertexIndex startIndex = 0;

      if (startFound == _coordinateMap.end()) {
        _coordinateMap.insert(std::make_pair(startCoord, _numVertices));
        startIndex = _numVertices;
        ++_numVertices;
      } else {
        startIndex = startFound->second;
      }

      VertexIndex endIndex = 0;

      if (endFound == _coordinateMap.end()) {
        _coordinateMap.insert(std::make_pair(endCoord, _numVertices));
        endIndex = _numVertices;
        ++_numVertices;
      } else {
        endIndex = endFound->second;
      }

      const std::pair<VertexIndex, VertexIndex> edge(startIndex, endIndex);

      _edgeMap.insert(
          std::make_pair(edge, std::make_pair(faceIndex, INVALID_INDEX)));

      // std::cerr << "face " << faceIndex << " edge " << edge.first << "-" <<
      // edge.second << "\n";
    }
  }
}

SurfaceGraph::SurfaceGraph(const PolyhedralSurface &surf)
    : _numVertices(0), _isValid(Validity::valid())
{
  const size_t numPatches = surf.numPatches();

  for (size_t p = 0; p != numPatches; ++p) { // for each polygon
    const FaceIndex idx = boost::add_vertex(_graph);
    BOOST_ASSERT(idx == p);
    (void)idx;
    const Polygon &polygon  = surf.patchN(p);
    const size_t   numRings = polygon.numRings();

    for (size_t r = 0; r != numRings; ++r) { // for each ring
      addRing(polygon.ringN(r), p);
    }
  }
}

SurfaceGraph::SurfaceGraph(const TriangulatedSurface &tin)
    : _numVertices(0), _isValid(Validity::valid())
{
  const size_t numPatches = tin.numPatches();

  for (size_t t = 0; t != numPatches; ++t) { // for each polygon
    const FaceIndex idx = boost::add_vertex(_graph);
    BOOST_ASSERT(idx == t);
    (void)idx;
    const Triangle &triangle = tin.patchN(t);
    addRing(triangle.toPolygon().exteriorRing(), t);
  }
}

auto
isConnected(const SurfaceGraph &graph) -> bool
{
  std::vector<SurfaceGraph::FaceIndex> component(
      boost::num_vertices(graph.faceGraph()));
  const size_t numComponents =
      boost::connected_components(graph.faceGraph(), component.data());
  return 1 == numComponents;
}

auto
isClosed(const SurfaceGraph &graph) -> bool
{
  const auto end = graph.edgeMap().end();

  for (auto e = graph.edgeMap().begin(); e != end; ++e) {
    if (e->second.second == SurfaceGraph::INVALID_INDEX) {
      return false;
    }
  }

  return true;
}

} // namespace SFCGAL::algorithm
