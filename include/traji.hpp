#ifndef TRAJI_HPP
#define TRAJI_HPP

#include <vector>
#include <utility>
#include <tuple>
#include <boost/geometry.hpp>

namespace traji {

typedef float TFloat;
typedef boost::geometry::model::point<TFloat, 2, boost::geometry::cs::cartesian> Point;
typedef boost::geometry::model::linestring<Point> LineString;

// forward declaration
class Path;
class Trajectory;

/// Representation of a position in the path (and possibly extended line)
/// When the position is on the path, 0 <= segment < #segments, 0 <= fraction < 1
/// When the position is on the extended last segment, segment = #segments - 1, fraction >= 1
/// When the position is on the extended first segment, segment = 0, fraction < 0
struct PathPosition
{
    std::size_t component; // index in multi-linestring
    std::size_t segment; // segment index in the line string
    TFloat fraction; // fraction of the point in the segment
    
    /// Convert the position to the distance to the beginning
    TFloat to_s(const Path &path);

    /// Convert the distance to the beginning to s
    static PathPosition from_s(const Path &path, TFloat s);
};

// TODO: this class should not be in public API
struct Segment
{
    /// Get the point that meets (p-p0) = fraction * (p1-p0)
    inline static Point interpolate(const Point &p0, const Point &p1, TFloat fraction)
    {
        auto x0 = p0.get<0>(), y0 = p0.get<1>();
        auto x1 = p0.get<0>(), y1 = p0.get<1>();
        return Point { x0 + (x1 - x0) * fraction, y0 + (y1 - y0) * fraction };
    }

    /// Calculate the signed distance from point to the segment **line**
    /// @return (distance to line, fraction of the foot point)
    /// The distance is positive if the point is at left hand side of the direction of line (p0 -> p1)
    static std::pair<TFloat, TFloat> sdistance(const Point &p0, const Point &p1, const TFloat l, const Point &rhs);

    /// Convert the sdistance to normal unsigned (squared) distance
    static inline TFloat distance2(const std::pair<TFloat, TFloat> &sdist, const TFloat l)
    {
        auto ds = sdist.first, d0 = sdist.second * l;
        if (sdist.second < 0)
            return ds * ds + d0 * d0;
        else if (sdist.second > 1)
            return ds * ds + (d0 - 1) * (d0 - 1);
        else
            return ds * ds;
    }
};

/// (immutable) non-parametric linestring
class Path
{
protected:
    LineString _line;
    std::vector<TFloat> _distance;

    void update_distance();

public:
    friend class PathPosition;

    inline TFloat length() const { return _distance.back(); }

    /// Get the point indicated by the distance from the beginning
    inline Point point_from(TFloat s) const
    {
        return point_at(PathPosition::from_s(*this, s));
    }
    Point point_at(const PathPosition &pos) const;

    /// Get the tagent at the point indicated by the distance from the beginning
    TFloat tangent_from(TFloat s) const
    {
        return tangent_at(PathPosition::from_s(*this, s));
    }
    TFloat tangent_at(const PathPosition &pos) const;

    /// Get the value interpolated by the distance from the beginning
    TFloat interpolate_from(const std::vector<TFloat> &values, TFloat s) const;

    /// Get the signed distance and foot point from a point to the path
    /// @return (distance to the path, projection point position)
    std::pair<TFloat, PathPosition> project(const Point &point) const;

    /// Access the underlying boost linestring
    inline const LineString& data() const { return _line; }

    /// Return the path in a form with equally space points
    /// @param smooth_radius If the radius is larger than zero, the corner of the polyline
    ///                      will be rounded by the radius when interpolating. If the radius is zero,
    ///                      the corner of the polyline will be reserved. Otherwise (by default), the
    ///                      corner of polyline is not reserved.
    Path respacing(TFloat resolution, float smooth_radius = -1) const;
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

    /// Return the trajectory represented by points with equal time intervals
    Trajectory reperiodize(TFloat resolution) const;
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
