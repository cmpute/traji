#ifndef TRAJI_HPP
#define TRAJI_HPP

#include <vector>
#include <utility>
#include <boost/geometry.hpp>

namespace traji {

typedef float TFloat;
typedef boost::geometry::model::point<TFloat, 2, boost::geometry::cs::cartesian> Point;
typedef boost::geometry::model::linestring<Point> LineString;

/// representation of position in the path
struct PathPosition
{
    std::size_t component;
    std::size_t segment;
    TFloat fraction;
};

/// (immutable) non-parametric linestring
class Path
{
protected:
    LineString _line;
    std::vector<TFloat> _distance;

public:
    inline TFloat length() const { return _distance.back(); }

    /// Get the point indicated by the distance from the beginning
    Point point_from(TFloat s) const;
    Point point_at(PathPosition &pos) const;

    /// Get the tagent at the point indicated by the distance from the beginning
    TFloat tangent_from(TFloat s) const;
    TFloat tangent_at(PathPosition &pos) const;

    /// Get the distance from the beginning. If the point is not on the line,
    /// then the closest point is used
    PathPosition project(const Point &point) const;

    /// Convert the position to the distance to the beginning
    TFloat to_s(const PathPosition &path);

    /// Get the value interpolated by the distance from the beginning
    TFloat interpolate_from(const std::vector<TFloat> &values, TFloat s) const;

    /// Get the signed distance from a point to the path
    std::pair<TFloat, TFloat> sdistance(const Point &point) const;

    inline const LineString& data() const { return _line; }
};

/// non-parametric linestring with time
class Trajectory : Path
{
protected:
    std::vector<float> _timestamp;

public:
    /// Get the point indicated by the time
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

TFloat distance(const Path &lhs, const Point &rhs);
std::vector<Point> intersection(const Path &lhs, const Path &rhs);
std::vector<std::pair<PathPosition, PathPosition>> arg_intersection(const Path &lhs, const Path &rhs);

Point cartesian_to_frenet(const Path &ref, const Point &point);
Point frenet_to_cartesian(const Path &ref, const Point &point);
Path cartesian_to_frenet(const Path &ref, const Path &path);
Path frenet_to_cartesian(const Path &ref, const Path &path);

}

#endif // TRAJI_HPP
