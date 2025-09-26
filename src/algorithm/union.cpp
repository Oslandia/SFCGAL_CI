// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/union.h"
#include "SFCGAL/algorithm/differencePrimitives.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/triangulate/triangulate2DZ.h"

#include <algorithm>
#include <cstdio>

#define DEBUG_OUT                                                              \
  if (0)                                                                       \
  std::cerr << __FILE__ << ":" << __LINE__ << " debug: "

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

struct EmptyPrimitive {};

enum PrimitiveType {
  PrimitivePoint   = 0,
  PrimitiveSegment = 1,
  PrimitiveSurface = 2,
  PrimitiveVolume  = 3,
  PrimitiveEmpty   = 4
};

template <int Dim>
struct Segment_d : detail::Segment_d<Dim>::Type {
  using PointType     = typename detail::Point_d<Dim>::Type;
  using SegmentType   = typename detail::Segment_d<Dim>::Type;
  using PointVector   = typename std::vector<PointType>;
  using SegmentVector = typename std::vector<SegmentType>;

  Segment_d(const SegmentType &s) : SegmentType(s) {}
  void
  splitAt(const PointType &p)
  {
    _split.push_back(p);
  }
  void
  remove(const SegmentType &s)
  {
    _split.push_back(s.source());
    _split.push_back(s.target());
    _remove.push_back(s);
  }

  template <class OutputIterator>
  [[nodiscard]] auto
  pieces(OutputIterator out) const -> OutputIterator
  {
    PointVector points(1, this->source());
    points.insert(points.end(), _split.begin(), _split.end());
    points.push_back(this->target());
    std::sort(points.begin() + 1, points.end() - 1,
              Nearer<PointType>(this->source()));

    for (auto p = points.begin(), q = p + 1; q != points.end(); ++p, q++) {
      if (*p != *q) {
        PointType const m = CGAL::midpoint(*p, *q);
        auto            r = _remove.begin();

        for (; r != _remove.end() && !r->has_on(m); ++r) {
          ;
        }

        if (r == _remove.end()) {
          *out++ = SegmentType(*p, *q);
        }
      }
    }

    return out;
  }

  [[nodiscard]] auto
  pieces() const -> SegmentVector
  {
    SegmentVector result;
    (void)this->pieces(std::back_inserter(result));
    return result;
  }

private:
  PointVector   _split;
  SegmentVector _remove;
};

template <int Dim>
struct Surface_d {};

template <>
struct Surface_d<3> : Triangle_3 {
  using PointVector   = std::vector<SFCGAL::Point_2>;
  using SegmentVector = std::vector<Segment_2>;
  using SurfaceVector = std::vector<PointVector>;

  Surface_d(const Triangle_3 &s) : Triangle_3(s), _plane(s.supporting_plane())
  {
    this->splitAt(s);
  }

  void
  splitAt(const SFCGAL::Point_3 &p)
  {
    //@note this is a degenerated segment, but works anyway
    _split.emplace_back(_plane.to_2d(p), _plane.to_2d(p));
  }

  void
  splitAt(const Segment_3 &s)
  {
    _split.emplace_back(_plane.to_2d(s.source()), _plane.to_2d(s.target()));
  }

  template <typename Point3Iterator>
  void
  splitAt(Point3Iterator begin, Point3Iterator end)
  { // polygon with unclosed ring
    if (begin == end) {
      return;
    }

    Point3Iterator s = begin;
    Point3Iterator t = s + 1;

    for (; t != end; ++t, ++s) {
      _split.push_back(Segment_2(_plane.to_2d(*s), _plane.to_2d(*t)));
    }

    _split.push_back(Segment_2(_plane.to_2d(*s), _plane.to_2d(*begin)));
  }

  void
  splitAt(const Triangle_3 &t)
  {
    const SFCGAL::Point_3 v[3] = {t.vertex(0), t.vertex(1), t.vertex(2)};
    this->splitAt(v, v + 3);
  }

  void
  splitAt(const std::vector<SFCGAL::Point_3> &p)
  { // polygon with unclosed ring
    this->splitAt(p.begin(), p.end());
  }

  template <typename Point3Iterator> // polygon with unclosed ring
  void
  remove(Point3Iterator begin, Point3Iterator end)
  {
    if (begin == end) {
      return;
    }

    this->splitAt(begin, end);
    PointVector v;

    for (Point3Iterator i = begin; i != end; ++i) {
      v.push_back(_plane.to_2d(*i));
    }

    _remove.push_back(v);
  }

  void
  remove(const std::vector<SFCGAL::Point_3> &p)
  {
    this->remove(p.begin(), p.end());
  }

  void
  remove(const Triangle_3 &t)
  {
    const SFCGAL::Point_3 v[3] = {t.vertex(0), t.vertex(1), t.vertex(2)};
    this->remove(v, v + 3);
  }

  auto
  pieces() -> std::vector<Triangle_3>
  {
    // we need to process the split lines because there may be several lines at
    // the same place and this won't play nice with the triangulation, the same
    // stands for the lines lying on the triangle edges
    //
    // after that we just check, for each triangle, if a point fall in a removed
    // part and remove it we can do that by pairwise union of all segments
    std::vector<Segment_2> filtered;
    {
      std::vector<Segment_d<2>> lines(_split.begin(), _split.begin() + 3);

      for (auto c = _split.begin() + 3; c != _split.end(); ++c) {
        Segment_d<2> current(*c);

        for (auto &line : lines) {
          CGAL::Object const inter = CGAL::intersection(line, current);
          const auto        *p     = CGAL::object_cast<Point_2>(&inter);
          const auto        *s     = CGAL::object_cast<Segment_2>(&inter);

          if (p != nullptr) {
            current.splitAt(*p);
            line.splitAt(*p);
          } else if (s != nullptr) {
            current.remove(*s);
            line.splitAt(s->source());
            line.splitAt(s->target());
          }
        }

        lines.push_back(current);
      }

      for (auto &line : lines) {
        (void)line.pieces(std::back_inserter(filtered));
      }
    }

    DEBUG_OUT << "triangulating " << filtered.size() << " lines\n";

    // now we want to triangulate
    TriangulatedSurface ts;
    {
      using Vertex_handle =
          triangulate::ConstraintDelaunayTriangulation::Vertex_handle;
      triangulate::ConstraintDelaunayTriangulation cdt;

      for (auto &f : filtered) {
        Vertex_handle const s = cdt.addVertex(f.source());
        Vertex_handle const t = cdt.addVertex(f.target());
        cdt.addConstraint(s, t);
      }

      cdt.getTriangles(ts);
    }

    // filter removed triangles
    std::vector<Triangle_3> res;

    for (auto &t : ts) {
      // define a point inside triangle
      const Point_2 a(t.vertex(0).toPoint_2());
      const Point_2 b(t.vertex(1).toPoint_2());
      const Point_2 c(t.vertex(2).toPoint_2());
      const Point_2 point(a + (Vector_2(a, b) + Vector_2(a, c)) / 3);

      // find if triangle is in a removed spot
      auto r = _remove.begin();

      for (; r != _remove.end() &&
             CGAL::bounded_side_2(r->begin(), r->end(), point, Kernel()) ==
                 CGAL::ON_UNBOUNDED_SIDE;
           ++r) {
        ;
      }

      if (r == _remove.end()) {
        res.emplace_back(_plane.to_3d(a), _plane.to_3d(b), _plane.to_3d(c));
      }
    }

    DEBUG_OUT << "generated " << res.size() << " triangles\n";

    return res;
  }

private:
  SFCGAL::Plane_3 _plane;
  SegmentVector   _split;
  SurfaceVector   _remove;
};

template <>
struct Surface_d<2> : Polygon_with_holes_2 {
  using PointVector   = std::vector<Point_2>;
  using SegmentVector = std::vector<Segment_2>;
  using SurfaceVector = std::vector<PointVector>;

  Surface_d(const Polygon_with_holes_2 &s) : Polygon_with_holes_2(s) {}

  void
  splitAt(const Segment_2 &s)
  {
    _split.push_back(s);
  }

  void
  addSplitsFrom(const Surface_d<2> &other)
  {
    _split.insert(_split.end(), other._split.begin(), other._split.end());
  }

  [[nodiscard]] auto
  pieces() const -> std::vector<Polygon_with_holes_2>
  {
    std::vector<Polygon_with_holes_2> res;
    fix_cgal_valid_polygon(*this, std::back_inserter(res));
    return res;
  }

private:
  SegmentVector _split;
};

// for debug prints
template <typename T>
auto
operator<<(std::ostream &out, std::set<T *> &obs) -> std::ostream &
{
  for (typename std::set<T *>::iterator h = obs.begin(); h != obs.end(); ++h) {
    out << *h << "\n";
  }

  return out;
}

// takes care of RAII of primitives

template <int Dim>
class Handle {
  struct ObservablePrimitive
      : boost::variant<typename detail::Point_d<Dim>::Type, Segment_d<Dim>,
                       Surface_d<Dim>, typename detail::Volume_d<Dim>::Type,
                       EmptyPrimitive> {
    template <class T>
    ObservablePrimitive(const T &p)
        : boost::variant<typename detail::Point_d<Dim>::Type, Segment_d<Dim>,
                         Surface_d<Dim>, typename detail::Volume_d<Dim>::Type,
                         EmptyPrimitive>(p)
    {
    }

    template <class T>
    auto
    as() -> T &
    {
      return boost::get<T &>(*this);
    }

    std::set<ObservablePrimitive **>
        _observers; // this is for ref counting and handle updating

  private:
    // non copyable
    ObservablePrimitive(const ObservablePrimitive &) = delete;
    auto
    operator=(const ObservablePrimitive &) -> ObservablePrimitive & = delete;
  };

public:
  Handle()
      : _p(new ObservablePrimitive *(new ObservablePrimitive(EmptyPrimitive())))
  {
    (*_p)->_observers.insert(_p);
    BOOST_ASSERT((*_p)->_observers.count(_p));
  }

  ~Handle()
  {
    BOOST_ASSERT((*_p)->_observers.count(_p));
    (*_p)->_observers.erase(_p);

    if ((*_p)->_observers.empty()) {
      delete (*_p);
    }

    delete _p;
  }

  Handle(const Handle &other) : _p(new ObservablePrimitive *(*other._p))
  {
    (*_p)->_observers.insert(_p);
    BOOST_ASSERT((*_p)->_observers.count(_p));
  }

  template <class PrimitiveType>
  explicit Handle(const PrimitiveType &primitive)
      : _p(new ObservablePrimitive *(new ObservablePrimitive(primitive)))
  {
    (*_p)->_observers.insert(_p);
    BOOST_ASSERT((*_p)->_observers.count(_p));
  }

  void
  swap(Handle &other) noexcept
  {
    (*_p)->_observers.erase(_p);
    (*other._p)->_observers.erase(other._p);
    std::swap(_p, other._p);
    (*_p)->_observers.insert(_p);
    (*other._p)->_observers.insert(other._p);
    BOOST_ASSERT((*_p)->_observers.count(_p));
    BOOST_ASSERT((*other._p)->_observers.count(other._p));
  }

  auto
  operator=(Handle other) -> Handle &
  {
    this->swap(other);
    BOOST_ASSERT((*_p)->_observers.count(_p));
    return *this;
  }

  auto
  operator*() const -> const ObservablePrimitive &
  {
    BOOST_ASSERT((*_p)->_observers.count(_p));
    return *(*_p);
  }
  auto
  operator*() -> ObservablePrimitive &
  {
    BOOST_ASSERT((*_p)->_observers.count(_p));
    return *(*_p);
  }
  auto
  operator->() const -> const ObservablePrimitive *
  {
    BOOST_ASSERT((*_p)->_observers.count(_p));
    return (*_p);
  }
  auto
  operator->() -> ObservablePrimitive *
  {
    BOOST_ASSERT((*_p)->_observers.count(_p));
    return (*_p);
  }

  auto
  asPoint() -> typename detail::Point_d<Dim>::Type &
  {
    BOOST_ASSERT((*_p)->_observers.count(_p));
    BOOST_ASSERT(which() == PrimitivePoint);
    return boost::get<typename detail::Point_d<Dim>::Type &>(*(*_p));
  }

  auto
  asSegment() -> Segment_d<Dim> &
  {
    BOOST_ASSERT((*_p)->_observers.count(_p));
    BOOST_ASSERT(which() == PrimitiveSegment);
    return boost::get<Segment_d<Dim> &>(*(*_p));
  }

  auto
  asSurface() -> Surface_d<Dim> &
  {
    BOOST_ASSERT((*_p)->_observers.count(_p));
    BOOST_ASSERT(which() == PrimitiveSurface);
    return boost::get<Surface_d<Dim> &>(*(*_p));
  }

  auto
  asVolume() -> typename detail::Volume_d<Dim>::Type &
  {
    BOOST_ASSERT((*_p)->_observers.count(_p));
    BOOST_ASSERT(which() == PrimitiveVolume);
    return boost::get<typename detail::Volume_d<Dim>::Type &>(*(*_p));
  }

  auto
  which() -> PrimitiveType
  {
    return PrimitiveType((*_p)->which());
  }

  auto
  empty() -> bool
  {
    return which() == PrimitiveEmpty;
  }

  // makes all handles observing a observe 'this' instead
  void
  registerObservers(Handle a)
  {
    if (*a._p == *_p) {
      return; // both aready observing the same primitive
    }

    ObservablePrimitive                *observed = *(a._p);
    std::vector<ObservablePrimitive **> observers(observed->_observers.begin(),
                                                  observed->_observers.end());

    for (auto h = observers.begin(); h != observers.end(); ++h) {
      *(*h) = *_p;
      (*(*h))->_observers.insert(*h);
    }

    BOOST_ASSERT(*a._p == *_p);

    delete observed; // we removed all observers
    BOOST_ASSERT((*_p)->_observers.count(_p));
#ifdef DEBUG

    for (typename std::vector<ObservablePrimitive **>::iterator h =
             observers.begin();
         h != observers.end(); ++h) {
      BOOST_ASSERT((*(*h))->_observers.count(*h));
    }

#endif
  }

private:
  ObservablePrimitive *
      *_p; // mutable because no non const cpy ctor en template classes
};

template <int Dim>
struct HandledBox {
  using Type = CGAL::Box_intersection_d::Box_with_handle_d<
      double, Dim, Handle<Dim>, CGAL::Box_intersection_d::ID_EXPLICIT>;
  using Vector = std::vector<Type>;
};

template <int Dim, class OutputIterator>
auto
compute_bboxes(const detail::GeometrySet<Dim> &gs, OutputIterator out)
    -> OutputIterator
{
  typename HandledBox<Dim>::Vector const bboxes;

  for (auto it = gs.points().begin(); it != gs.points().end(); ++it) {
    *out++ = typename HandledBox<Dim>::Type(it->primitive().bbox(),
                                            Handle<Dim>(it->primitive()));
  }

  for (auto it = gs.segments().begin(); it != gs.segments().end(); ++it) {
    *out++ = typename HandledBox<Dim>::Type(it->primitive().bbox(),
                                            Handle<Dim>(it->primitive()));
  }

  for (auto it = gs.surfaces().begin(); it != gs.surfaces().end(); ++it) {
    DEBUG_OUT << "Adding surface " << it->primitive() << "\n";
    DEBUG_OUT << "       surface box " << it->primitive().bbox() << "\n";
    *out++ = typename HandledBox<Dim>::Type(it->primitive().bbox(),
                                            Handle<Dim>(it->primitive()));
  }

  for (auto it = gs.volumes().begin(); it != gs.volumes().end(); ++it) {
    *out++ = typename HandledBox<Dim>::Type(
        compute_solid_bbox(it->primitive(), detail::dim_t<Dim>()),
        Handle<Dim>(it->primitive()));
  }

  return out;
}

template <class Handle>
void
union_point_point(Handle a, Handle b)
{
  if (a.asPoint() == b.asPoint()) {
    a.registerObservers(b);
  }
}

template <class Handle>
void
union_point_segment(Handle a, Handle b)
{
  if (b.asSegment().has_on(a.asPoint())) {
    b.asSegment().splitAt(a.asPoint());
    b.registerObservers(a);
  }
}

void
union_point_surface(Handle<2> a, Handle<2> b)
{
  if (do_intersect(a.asPoint(), b.asSurface())) {
    b.registerObservers(a);
  }
}

void
union_point_surface(Handle<3> a, Handle<3> b)
{
  if (b.asSurface().has_on(a.asPoint())) {
    b.registerObservers(a);
  }
}

void
union_point_volume(const Handle<2> & /*unused*/, const Handle<2> & /*unused*/)
{
  BOOST_ASSERT(false); // there shouldn't be any volume in 2D
}

void
union_point_volume(Handle<3> a, Handle<3> b)
{
  //@todo put is in poly in a struct derived from MarkedPolyhedron to avoid
  // rebuilding point inside every time
  CGAL::Side_of_triangle_mesh<MarkedPolyhedron, Kernel> const is_in_poly(
      b.asVolume());

  if (CGAL::ON_UNBOUNDED_SIDE != is_in_poly(a.asPoint())) {
    b.registerObservers(a);
  }
}

template <int Dim>
void
union_segment_segment(Handle<Dim> a, Handle<Dim> b)
{
  using SegmentType = typename detail::TypeForDimension<Dim>::Segment;
  using PointType   = typename detail::TypeForDimension<Dim>::Point;

  CGAL::Object const inter = CGAL::intersection(a.asSegment(), b.asSegment());
  const auto        *p     = CGAL::object_cast<PointType>(&inter);
  const auto        *s     = CGAL::object_cast<SegmentType>(&inter);

  if (p) {
    b.asSegment().splitAt(*p);
    a.asSegment().splitAt(*p);
  } else if (s) {
    b.asSegment().remove(*s);
    a.asSegment().splitAt(s->source());
    a.asSegment().splitAt(s->target());
  }
}

void
union_segment_segment(const Handle<2> &a, const Handle<2> &b)
{
  union_segment_segment<2>(a, b);
}

void
union_segment_segment(const Handle<3> &a, const Handle<3> &b)
{
  union_segment_segment<3>(a, b);
}

void
union_segment_surface(Handle<2> a, Handle<2> b)
{
  std::vector<Polygon_2> rings(1, b.asSurface().outer_boundary());
  rings.insert(rings.end(), b.asSurface().holes_begin(),
               b.asSurface().holes_end());

  std::vector<Point_2> points(1, a.asSegment().source());

  for (auto &ring : rings) {
    for (auto target = ring.vertices_begin(); target != ring.vertices_end();
         ++target) {
      const Segment_2 sc(target == ring.vertices_begin()
                             ? *(ring.vertices_end() - 1)
                             : *(target - 1),
                         *target);

      CGAL::Object const inter = CGAL::intersection(a.asSegment(), sc);
      const auto        *p     = CGAL::object_cast<Point_2>(&inter);
      const auto        *s     = CGAL::object_cast<Segment_2>(&inter);

      if (p != nullptr) {
        points.push_back(*p);
      } else if (s != nullptr) {
        a.asSegment().remove(*s);
      }
    }
  }

  points.push_back(a.asSegment().target());

  // order point according to the distance from source
  std::sort(points.begin() + 1, points.end() - 1, Nearer<Point_2>(points[0]));

  // cut segment with pieces that have length and wich midpoint is inside
  // polygon
  for (auto p = points.begin(), q = p + 1; q != points.end(); ++p, ++q) {
    if (*p != *q && do_intersect(CGAL::midpoint(*p, *q), b.asSurface())) {
      const Segment_2 s(*p, *q);
      a.asSegment().remove(s);
      b.asSurface().splitAt(s);
    }
  }
}

void
union_segment_surface(Handle<3> a, Handle<3> b)
{
  CGAL::Object const inter = CGAL::intersection(a.asSegment(), b.asSurface());
  const auto        *s     = CGAL::object_cast<Segment_3>(&inter);

  if (s != nullptr) {
    a.asSegment().remove(*s);
    b.asSurface().splitAt(*s);
  }
}

void
union_segment_volume(const Handle<2> & /*unused*/, const Handle<2> & /*unused*/)
{
  BOOST_ASSERT(false); // there shouldn't be any volume in 2D
}

void
union_segment_volume(Handle<3> a, Handle<3> b)
{
  const Segment_3        &segment    = a.asSegment();
  const MarkedPolyhedron &polyhedron = b.asVolume();

  std::vector<FaceBbox>     bboxes(polyhedron.facets_begin(),
                                   polyhedron.facets_end());
  std::vector<FaceBboxBase> bbox(
      1, FaceBboxBase(segment.bbox(),
                      polyhedron.facets_begin()->facet_begin())); // nevermind
                                                                  // the facet
                                                                  // handle,
                                                                  // it's not
                                                                  // used anyway
  FaceSegmentCollide::CollisionVector collisions;
  FaceSegmentCollide const            cb(collisions);
  CGAL::box_intersection_d(bbox.begin(), bbox.end(), bboxes.begin(),
                           bboxes.end(), cb);

  CGAL::Side_of_triangle_mesh<MarkedPolyhedron, Kernel> const is_in_poly(
      polyhedron);

  if (collisions.empty()) {
    // completely in or out, we just test one point

    if (CGAL::ON_UNBOUNDED_SIDE != is_in_poly(segment.source())) {
      a.asSegment().remove(a.asSegment());
    }
  } else {
    std::vector<Triangle_3> triangles;
    collidingTriangles(collisions, std::back_inserter(triangles));

    // first step, substract faces
    for (auto &triangle : triangles) {
      Handle<3> const h(triangle);
      union_segment_surface(a, h);
    }

    // second step, for each segment, add intersection points and test each
    // middle point to know if it's in or out
    std::vector<Point_3> points;

    for (auto &triangle : triangles) {
      CGAL::Object const inter = CGAL::intersection(segment, triangle);
      const auto        *p     = CGAL::object_cast<Point_3>(&inter);

      if (p != nullptr) {
        points.push_back(*p);
      }
    }

    if (static_cast<unsigned int>(!points.empty()) != 0U) {
      std::sort(points.begin(), points.end(),
                Nearer<Point_3>(segment.source()));

      // mark segments pieces that have length and wich midpoint is inside
      // polyhedron
      for (auto p = points.begin(), q = p + 1; q != points.end(); ++p, ++q) {
        if (*p != *q &&
            CGAL::ON_UNBOUNDED_SIDE != is_in_poly(CGAL::midpoint(*p, *q))) {
          a.asSegment().remove(Segment_3(*p, *q));
        }
      }
    }
  }
}

void
union_surface_surface(Handle<2> a, Handle<2> b)
{
  Polygon_with_holes_2 res;

  if (CGAL::join(fix_sfs_valid_polygon(a.asSurface()),
                 fix_sfs_valid_polygon(b.asSurface()), res)) {
    DEBUG_OUT << "merged " << a.asSurface() << " and " << b.asSurface() << "\n";
    Handle<2> h(res);
    h.asSurface().addSplitsFrom(a.asSurface());
    h.asSurface().addSplitsFrom(b.asSurface());
    h.registerObservers(a);
    h.registerObservers(b);
  }
}

void
union_surface_surface(Handle<3> a, Handle<3> b)
{
  CGAL::Object const inter = intersection(a.asSurface(), b.asSurface());
  const auto        *p     = CGAL::object_cast<Point_3>(&inter);
  const auto        *s     = CGAL::object_cast<Segment_3>(&inter);
  const auto        *t     = CGAL::object_cast<Triangle_3>(&inter);
  const auto        *v     = CGAL::object_cast<std::vector<Point_3>>(&inter);

  if (p != nullptr) {
    a.asSurface().splitAt(*p);
    b.asSurface().splitAt(*p);
  } else if (s != nullptr) {
    a.asSurface().splitAt(*s);
    b.asSurface().splitAt(*s);
  } else if (t != nullptr) {
    a.asSurface().splitAt(*t);
    b.asSurface().remove(*t);
  } else if (v != nullptr) {
    a.asSurface().splitAt(*v);
    b.asSurface().remove(*v);
  }
}

void
union_surface_volume(const Handle<2> & /*unused*/, const Handle<2> & /*unused*/)
{
  BOOST_ASSERT(false); // there shouldn't be any volume in 2D
}

void
union_surface_volume(Handle<3> a, Handle<3> b)
{
  detail::GeometrySet<3> res;
  intersectionSolidTriangle(b.asVolume(), a.asSurface(), res);

  for (auto &it : res.surfaces()) {
    a.asSurface().remove(it.primitive());
  }
}

void
union_volume_volume(const Handle<2> & /*unused*/, const Handle<2> & /*unused*/)
{
  BOOST_ASSERT(false); // there shouldn't be any volume in 2D
}

void
union_volume_volume(Handle<3> a, Handle<3> b)
{
  auto &p = const_cast<MarkedPolyhedron &>(a.asVolume());
  auto &q = const_cast<MarkedPolyhedron &>(b.asVolume());

  // volumes must at least share a face, if they share only a point, this will
  // cause an invalid geometry, if they only share an egde it will cause the
  // CGAL algo to throw
  //
  // for the moment we use hte intersection algo and test the result
  detail::GeometrySet<3> inter;
  intersection(detail::GeometrySet<3>(a.asVolume()),
               detail::GeometrySet<3>(b.asVolume()), inter);

  if ((static_cast<unsigned int>(!inter.volumes().empty()) != 0U) ||
      (static_cast<unsigned int>(!inter.surfaces().empty()) != 0U)) {

    MarkedPolyhedron output;
    bool const       res =
        CGAL::Polygon_mesh_processing::corefine_and_compute_union(p, q, output);

    if (res && std::next(vertices(output).first) != vertices(output).second) {
      Handle<3> h(output);
      // @todo check that the volume is valid (connection on one point isn't)
      h.registerObservers(a);
      h.registerObservers(b);
    }
  }
}

template <int Dim>
struct UnionOnBoxCollision {
  void
  operator()(typename HandledBox<Dim>::Type &a,
             typename HandledBox<Dim>::Type &b)
  {
    DEBUG_OUT << "collision of boxes\n";

    switch (a.handle().which()) {
    case PrimitivePoint:
      switch (b.handle().which()) {
      case PrimitivePoint:
        union_point_point(a.handle(), b.handle());
        break;

      case PrimitiveSegment:
        union_point_segment(a.handle(), b.handle());
        break;

      case PrimitiveSurface:
        union_point_surface(a.handle(), b.handle());
        break;

      case PrimitiveVolume:
        union_point_volume(a.handle(), b.handle());
        break;

      case PrimitiveEmpty:
        break;
      }

      break;

    case PrimitiveSegment:
      switch (b.handle().which()) {
      case PrimitivePoint:
        union_point_segment(b.handle(), a.handle());
        break;

      case PrimitiveSegment:
        union_segment_segment(a.handle(), b.handle());
        break;

      case PrimitiveSurface:
        union_segment_surface(a.handle(), b.handle());
        break;

      case PrimitiveVolume:
        union_segment_volume(a.handle(), b.handle());
        break;

      case PrimitiveEmpty:
        break;
      }

      break;

    case PrimitiveSurface:
      switch (b.handle().which()) {
      case PrimitivePoint:
        union_point_surface(b.handle(), a.handle());
        break;

      case PrimitiveSegment:
        union_segment_surface(b.handle(), a.handle());
        break;

      case PrimitiveSurface:
        union_surface_surface(a.handle(), b.handle());
        break;

      case PrimitiveVolume:
        union_surface_volume(a.handle(), b.handle());
        break;

      case PrimitiveEmpty:
        break;
      }

      break;

    case PrimitiveVolume:
      switch (b.handle().which()) {
      case PrimitivePoint:
        union_point_volume(b.handle(), a.handle());
        break;

      case PrimitiveSegment:
        union_segment_volume(b.handle(), a.handle());
        break;

      case PrimitiveSurface:
        union_surface_volume(b.handle(), a.handle());
        break;

      case PrimitiveVolume:
        union_volume_volume(a.handle(), b.handle());
        break;

      case PrimitiveEmpty:
        break;
      }

      break;

    case PrimitiveEmpty:
      break;
    }
  }
};

template <int Dim>
void
collectPrimitives(const typename HandledBox<Dim>::Vector &boxes,
                  detail::GeometrySet<Dim>               &output)
{
  Handle<Dim> empty;

  for (auto bit = boxes.begin(); bit != boxes.end(); ++bit) {
    switch (bit->handle().which()) {
    case PrimitivePoint:
      output.addPrimitive(bit->handle().asPoint());
      empty.registerObservers(bit->handle());
      break;

    case PrimitiveSegment: {
      typename std::vector<typename detail::Segment_d<Dim>::Type> pieces(
          bit->handle().asSegment().pieces());
      output.addSegments(pieces.begin(), pieces.end());
      empty.registerObservers(bit->handle());
    } break;

    case PrimitiveSurface: {
      typename std::vector<typename detail::Surface_d<Dim>::Type> pieces(
          bit->handle().asSurface().pieces());
      output.addSurfaces(pieces.begin(), pieces.end());
      empty.registerObservers(bit->handle());
    } break;

    case PrimitiveVolume:
      output.addPrimitive(bit->handle().asVolume());
      empty.registerObservers(bit->handle());
      break;

    case PrimitiveEmpty:
      break;
    }
  }
}

auto
union_(const Geometry &ga, const Geometry &gb, NoValidityCheck /*unused*/)
    -> std::unique_ptr<Geometry>
{
  HandledBox<2>::Vector boxes;
  compute_bboxes(detail::GeometrySet<2>(ga), std::back_inserter(boxes));
  const unsigned numBoxA = boxes.size();
  compute_bboxes(detail::GeometrySet<2>(gb), std::back_inserter(boxes));

  CGAL::box_intersection_d(boxes.begin(), boxes.begin() + numBoxA,
                           boxes.begin() + numBoxA, boxes.end(),
                           UnionOnBoxCollision<2>());

  detail::GeometrySet<2> output;
  collectPrimitives(boxes, output);
  return output.recompose();
}

auto
union_(const Geometry &ga, const Geometry &gb) -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(ga);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(gb);
  std::unique_ptr<Geometry> result(union_(ga, gb, NoValidityCheck()));
  return result;
}

/// @} end of private section

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

/**
 * @brief Compute 3D union of two geometries without validity check
 * @param ga First geometry
 * @param gb Second geometry
 * @param NoValidityCheck Tag to skip validity checks
 * @return The 3D union of ga and gb
 */
auto
union3D(const Geometry &ga, const Geometry &gb, NoValidityCheck /*unused*/)
    -> std::unique_ptr<Geometry>
{
  HandledBox<3>::Vector boxes;
  compute_bboxes(detail::GeometrySet<3>(ga), std::back_inserter(boxes));
  const unsigned numBoxA = boxes.size();
  compute_bboxes(detail::GeometrySet<3>(gb), std::back_inserter(boxes));

  CGAL::box_intersection_d(boxes.begin(), boxes.begin() + numBoxA,
                           boxes.begin() + numBoxA, boxes.end(),
                           UnionOnBoxCollision<3>());

  detail::GeometrySet<3> output;
  collectPrimitives(boxes, output);
  return output.recompose();
}

/**
 * @brief Compute 3D union of two geometries
 * @param ga First geometry
 * @param gb Second geometry
 * @return The 3D union of ga and gb
 */
auto
union3D(const Geometry &ga, const Geometry &gb) -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(ga);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(gb);
  std::unique_ptr<Geometry> result(union3D(ga, gb, NoValidityCheck()));
  return result;
}

/// @private
void
handleLeakTest()
{
  Handle<2> const h0(Point_2(0, 0));
  Handle<2>       h1(Point_2(1, 1));
  Handle<2>       empty;
  empty.registerObservers(empty);
  empty.registerObservers(h0);
  h1.registerObservers(h0);
}
} // namespace SFCGAL::algorithm
