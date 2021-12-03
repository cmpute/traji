#ifndef TRAJI_HPP
#define TRAJI_HPP

#include <vector>
#include <utility>
#include <tuple>
#include <boost/geometry.hpp>
#include <Eigen/Core>

namespace traji {

typedef float TFloat;
typedef Eigen::Matrix<TFloat, 2, 1> Vector2;
typedef Eigen::Matrix<TFloat, 3, 1> Vector3;
typedef Eigen::Matrix<TFloat, 3, 3> Matrix3;
typedef Eigen::Matrix<TFloat, 6, 1> Vector6;
typedef Eigen::Matrix<TFloat, -1, 1> VectorX;
typedef Eigen::Array<TFloat, -1, 1> ArrayX;
typedef boost::geometry::model::point<TFloat, 2, boost::geometry::cs::cartesian> Point;
typedef boost::geometry::model::linestring<Point> LineString;

#ifndef M_PI
constexpr TFloat pi = 3.14159265358979323846;
constexpr TFloat pi2 = pi/2;
#else
constexpr TFloat pi = M_PI;
constexpr TFloat pi2 = M_PI_2;
#endif

// forward declaration
class Path;
class Trajectory;
class HeteroSegment;
class HeteroPath;
class QuinticPolyTrajectory;

enum class SegmentType
{
    Line, // no params
    Arc, // param: x, y of the arc center 
    QuadraticBezier, // param: x, y of the control point
    CubicBezier, // param: x1, y1, x2, y2 of the two control points

    // For polynomials, the latent parameter range need to be normalized to [0, 1]
    Polynomial, // params: polynomial coeffs (high to low)
};

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
    TFloat to_s(const HeteroPath &path);

    /// Convert the distance to the beginning to s
    static PathPosition from_s(const Path &path, TFloat s);
    static PathPosition from_s(const HeteroPath &path, TFloat s);
    static std::vector<PathPosition> from_s(const Path &path, const std::vector<TFloat> &s_list);

    /// Convert the position to the distance to the timestamp
    TFloat to_t(const Trajectory &traj);
    static PathPosition from_t(const Trajectory &traj, TFloat t);
};

/// (immutable) non-parametric linestring
// TODO: implement move semantics
class Path
{
protected:
    LineString _line; // the geometry of the line string
    std::vector<TFloat> _distance; // the distance of each vertices to the first vertex along the line. First value is always zero

    /// Update _distance values
    void update_distance();

    /// Calculate feasible smooth radius given the input as max limit
    std::vector<TFloat> calc_feasible_radius(TFloat smooth_radius) const;

private:
    /// respacing with smooth_radius = 0
    Path respacing0(TFloat resolution) const;

public:
    friend class PathPosition;
    friend class QuinticPolyTrajectory;

    Path() {}
    template<typename Iterator>
    Path(Iterator begin, Iterator end) : _line(begin, end) { update_distance(); }
    Path(std::initializer_list<Point> l) : _line(l) { update_distance(); }

    inline std::size_t size() const { return _line.size(); }
    inline TFloat length() const { return _distance.back(); }
    inline std::vector<TFloat> segment_lengths() const
    {
        std::vector<TFloat> lengths(_line.size() - 1);
        for (int i = 1; i < _line.size(); i++)
            lengths[i-1] = _distance[i] - _distance[i-1];
        return lengths;
    }

    /// Access the underlying boost linestring
    inline LineString& data() { return _line; }
    inline const LineString& data() const { return _line; }

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

    /// Return the path in a form with equally space points (interpolating over S)
    /// The starting and end point of the original line string is guaranteed to be included and not smoothed
    /// @param smooth_radius The distance from original vertex to the new smoothed vertex.
    Path respacing(TFloat resolution, TFloat smooth_radius = 0) const;

    /// Return the path with each segment's length less than the given resolution
    /// This method will retain the original corner points, unlike respacing function
    Path densify(TFloat resolution) const;

    /// Return the path with rounded corners. This method doesn't change the density of
    /// the points on the line segments
    HeteroPath smooth(TFloat smooth_radius, SegmentType segtype = SegmentType::Arc) const;

    /// Return a path represented by a list of distance to start. 
    /// The best performance is achieved when s_list is sorted ascendingly.
    Path resample(const std::vector<TFloat> &s_list) const;
};

/// non-parametric linestring with time
class Trajectory : public Path
{
protected:
    std::vector<TFloat> _timestamps; // the timestamp over each point

public:
    friend class PathPosition;
    friend class QuinticPolyTrajectory;

    Trajectory() {}
    Trajectory(const Path& path, const std::vector<TFloat> &timestamps)
        : Path(path), _timestamps(timestamps) {}
    template<typename Iterator, typename TIterator>
    Trajectory(Iterator begin, Iterator end, TIterator t_begin, TIterator t_end)
        : Path(begin, end), _timestamps(t_begin, t_end) { update_distance(); }

    const std::vector<TFloat>& timestamps() const { return _timestamps; }

    /// Get the point indicated by the time
    Point point_at(TFloat t) const;
    TFloat tangent_at(TFloat t) const;

    /// Return the trajectory represented by points with equal time intervals
    Trajectory reperiodize(TFloat interval) const;
};


// ==================================== Parametric paths ====================================

class QuinticPolyTrajectory
{
protected:
    Vector6 _x_coeffs, _y_coeffs; // high-order to low-order
    TFloat _T; // the trajectory starts at t=0 and ends at t=T

public:
    QuinticPolyTrajectory(TFloat T, const Vector6 &x_coeffs, const Vector6 &y_coeffs)
        : _x_coeffs(x_coeffs), _y_coeffs(y_coeffs), _T(T) {}

    /// Solve optimal trajectory based on x and y states (including 1st and 2nd state derivatives)
    QuinticPolyTrajectory(TFloat T, const Vector3 &x0, const Vector3 &xT,
                          const Vector3 &y0, const Vector3 &yT);

    const Vector6& x_coeffs() const { return _x_coeffs; }
    const Vector6& y_coeffs() const { return _y_coeffs; }
    TFloat T() const { return _T; }

    Point point_at(TFloat t) const;
    TFloat tangent_at(TFloat t) const;
    Vector2 velocity_at(TFloat t) const;

    Trajectory rasterize(TFloat resolution) const;
    Trajectory periodize(TFloat interval) const;
};

// TODO: spline interpolated paths are also parametric paths

// ==================================== Hybrid paths ====================================

struct HeteroSegment
{
    SegmentType type;
    Point start, end;

    /// Params for the segments, see SegmentType for details
    std::vector<TFloat> params;

    Point point_at(TFloat fraction);
    TFloat tangent_at(TFloat fraction);
};

class HeteroPath
{
protected:
    std::vector<HeteroSegment> _segments;
    std::vector<TFloat> _distance; // distance to the end of the segment

public:
    HeteroPath() {}
    HeteroPath(const Path& path);

    friend class Path;

    /// Convert to Path by rasterizing each segments
    Path rasterize(TFloat resolution);
    /// Convert to Path by only rasterizing curves (excl. line segment)
    Path rasterize_curve(TFloat resolution);

    Path respacing(TFloat resolution) const;
    HeteroPath densify(TFloat resolution) const;
};

// ==================================== frenet paths ====================================

/// In frenet coordinate, Point::x represents s along the path (forward = positive), Point::y
/// represents l along the normal direction of the path (left = positive)
namespace frenet
{
    Point from_cartesian(const Path &ref, const Point &point);
    Point to_cartesian(const Path &ref, const Point &point);

    // TODO: add conversion for dynamic states from/to frenet
    // REF: https://blog.csdn.net/davidhopper/article/details/79162385#t6

    Path from_cartesian(const Path &ref, const Path &path);
    Path to_cartesian(const Path &ref, const Path &path);

    inline Trajectory from_cartesian(const Path &ref, const Trajectory &traj)
    {
        return Trajectory(from_cartesian(ref, traj), traj.timestamps());
    }
    inline Trajectory to_cartesian(const Path &ref, const Trajectory &traj)
    {
        return Trajectory(to_cartesian(ref, traj), traj.timestamps());
    }
}

// ==================================== binary functions ====================================

inline TFloat distance(const Point &lhs, const Point &rhs) { return boost::geometry::distance(lhs, rhs); }
inline TFloat distance(const Path &lhs, const Point &rhs) { return lhs.project(rhs).first; }
inline PathPosition arg_distance(const Path &lhs, const Point &rhs) { return lhs.project(rhs).second; }

std::vector<Point> intersection(const Path &lhs, const Path &rhs);
std::vector<std::pair<PathPosition, PathPosition>> arg_intersection(const Path &lhs, const Path &rhs);

// ==================================== standard functions ====================================

inline bool operator==(const PathPosition& lhs, const PathPosition& rhs)
{
    return lhs.component == rhs.component && lhs.segment == rhs.segment && lhs.fraction == rhs.fraction;
}
inline bool operator!=(const PathPosition& lhs, const PathPosition& rhs) { return !(lhs == rhs); }

bool operator==(const Path& lhs, const Path& rhs);
inline bool operator!=(const Path& lhs, const Path& rhs) { return !(lhs == rhs); }

} // namespace traji

namespace boost { namespace geometry { namespace model
{
    
inline bool operator==(const traji::Point& lhs, const traji::Point& rhs)
{
    return lhs.get<0>() == rhs.get<0>() && lhs.get<1>() == rhs.get<1>();
}
inline bool operator!=(const traji::Point& lhs, const traji::Point& rhs) { return !(lhs == rhs); }

}}} // namespace boost.geometry.model

namespace std
{

std::string to_string(const traji::Point &value);
std::string to_string(const traji::Path &value);
std::string to_string(const traji::PathPosition &value);
std::string to_string(const traji::Trajectory &value);
std::string to_string(const traji::QuinticPolyTrajectory &value);

}

#endif // TRAJI_HPP
