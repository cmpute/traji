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

// ==================================== Non-parametric paths ====================================

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

/// (immutable) non-parametric linestring
class Path
{
protected:
    LineString _line; // the geometry of the line string
    std::vector<TFloat> _distance; // the distance of each vertices to the first vertex along the line

    /// Update _distance values
    void update_distance();

public:
    friend class PathPosition;

    Path() {}
    template<typename Iterator>
    Path(Iterator begin, Iterator end) : _line(begin, end) { update_distance(); }
    Path(std::initializer_list<Point> l) : _line(l) { update_distance(); }

    inline std::size_t size() const { return _line.size(); }
    inline TFloat length() const { return _distance.back(); }

    Point& operator[](std::size_t idx) { return _line[idx]; }
    const Point& operator[](std::size_t idx) const { return _line[idx]; }

    /// Get the point indicated by the distance from the beginning
    inline Point point_from(TFloat s) const
    {
        return point_at(PathPosition::from_s(*this, s));
    }
    Point point_at(const PathPosition &pos) const;

    /// Get the tagent at the point indicated by the distance from the beginning
    inline TFloat tangent_from(TFloat s) const
    {
        return tangent_at(PathPosition::from_s(*this, s));
    }
    TFloat tangent_at(const PathPosition &pos) const;

    /// Get the value interpolated by the distance from the beginning
    inline TFloat interpolate_from(const std::vector<TFloat> &values, TFloat s) const
    {
        return interpolate_at(values, PathPosition::from_s(*this, s));
    }
    TFloat interpolate_at(const std::vector<TFloat> &values, const PathPosition &pos) const;

    /// Get the signed distance and foot point from a point to the path
    /// @return (distance to the path, projection point position)
    std::pair<TFloat, PathPosition> project(const Point &point) const;

    /// Access the underlying boost linestring
    inline LineString& data() { return _line; }
    inline const LineString& data() const { return _line; }

    /// Return the path in a form with equally space points (interpolating over S)
    /// @param smooth_radius If the radius is larger than zero, the corner of the polyline
    ///                      will be rounded by a bezier curve with the given radius when interpolating.
    Path respacing(TFloat resolution, float smooth_radius = 0) const;

    /// Return the path with each segment's length less than the given resolution
    /// This method will retain the original corner points, unlike respacing function
    Path densify(TFloat resolution) const;
};

/// non-parametric linestring with time
class Trajectory : Path
{
protected:
    std::vector<float> _timestamp; // the timestamp over each point

public:
    /// Get the point indicated by the time
    Point point_at(TFloat t) const;
    TFloat tangent_at(TFloat t) const;

    /// Return the trajectory represented by points with equal time intervals
    Trajectory reperiodize(TFloat resolution) const;
};


// ==================================== Parametric paths ====================================

class QuinticPolynomialTrajectory
{
public:
    Point point_from(TFloat s) const;
    TFloat tangent_from(TFloat s) const;
    

    TFloat project(Point point) const;
    Trajectory rollout() const;
};

// ==================================== Hybrid paths ====================================

// TODO: support linestring with heterogeneous segment types (called HeteroPath?)
// this could be the better implementation for rounded corner
enum class SegmentType
{
    Line, // no params
    Arc, // param: radius
    QuadraticBezier, // param: x, y of the control point
    CubicBezier, // param: x1, y1, x2, y2 of the two control points

    // For polynomials, the argument range need to be normalized to [0, 1]
    QuarticPoly, // params: polynomial coeffs
    QuinticPoly // params: polynomial coeffs
};

struct HeteroSegment
{
    SegmentType type;
    Point start, end;
    std::vector<TFloat> params;
};

// ==================================== binary functions ====================================

inline TFloat distance(const Point &lhs, const Point &rhs) { return boost::geometry::distance(lhs, rhs); }
TFloat distance(const Path &lhs, const Point &rhs);

std::vector<Point> intersection(const Path &lhs, const Path &rhs);
std::vector<std::pair<PathPosition, PathPosition>> arg_intersection(const Path &lhs, const Path &rhs);

Point cartesian_to_frenet(const Path &ref, const Point &point);
Point frenet_to_cartesian(const Path &ref, const Point &point);
Path cartesian_to_frenet(const Path &ref, const Path &path);
Path frenet_to_cartesian(const Path &ref, const Path &path);

} // namespace traji


// ==================================== extended std functions ====================================

namespace std
{

std::string to_string(const traji::Point &value);
std::string to_string(const traji::Path &value);
std::string to_string(const traji::PathPosition &value);

}

#endif // TRAJI_HPP
