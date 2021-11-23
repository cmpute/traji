#ifndef TRAJI_HPP
#define TRAJI_HPP

#include <boost/geometry.hpp>

namespace traji {

typedef float TFloat;
typedef boost::geometry::model::point<TFloat, 2, boost::geometry::cs::cartesian> Point;
typedef boost::geometry::model::linestring<Point> LineString;

/// non-parametric linestring
class Path
{
protected:
    LineString _line;
public:
    Point point_from(TFloat s) const;
    TFloat tangent_from(TFloat s) const;
    TFloat project(Point point) const;
    // TODO: support interpolation of any associated field?
};

/// non-parametric linestring with time
class Trajectory : Path
{
protected:
    std::vector<float> timestamp;
public:
    Point point_at(TFloat t) const;
    TFloat tangent_at(TFloat t) const;
};

// parametric path / trajectory
class QuintPoly
{
public:
    Point point_from(TFloat s) const;
    TFloat tangent_from(TFloat s) const;
    TFloat project(Point point) const;
    Trajectory rollout() const;
};

Point cartesian_to_frenet(Path ref, Point point);
Point frenet_to_cartesian(Path ref, Point point);
Path cartesian_to_frenet(Path ref, Path path);
Path frenet_to_cartesian(Path ref, Path path);

}

#endif // TRAJI_HPP
