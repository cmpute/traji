#ifndef TRAJI_HPP
#define TRAJI_HPP

// TODO: some big changes
// 1. fix TRel to float, while always use double for s0 and t0, and store the difference from s0 / t0 instead.

#include <vector>
#include <utility>
#include <tuple>
#include <boost/geometry.hpp>
#include <Eigen/Core>

/* forward declarations */
namespace traji {

/// @brief Float type used for relative coordinates and values
typedef float TRel;

/// @brief Float type used for absolute coordinates
typedef double TAbs;

typedef Eigen::Matrix<TRel, 2, 1> Vector2;
typedef Eigen::Matrix<TRel, 3, 1> Vector3;
typedef Eigen::Matrix<TRel, 2, 2> Matrix2;
typedef Eigen::Matrix<TRel, 3, 3> Matrix3;
typedef Eigen::Matrix<TRel, 6, 1> Vector6;
typedef Eigen::Matrix<TRel, -1, 1> VectorX;
typedef Eigen::Array<TRel, -1, 1> ArrayX;

#ifndef M_PI
constexpr TRel pi = 3.14159265358979323846;
constexpr TRel pi2 = pi/2;
#else
constexpr TRel pi = M_PI;
constexpr TRel pi2 = M_PI_2;
#endif

typedef boost::geometry::model::point<TRel, 2, boost::geometry::cs::cartesian> Point;
typedef boost::geometry::model::linestring<Point> LineString;

class Path;
class Trajectory;
class HeteroSegment;
class HeteroPath;
class QuinticPolyTrajectory;
class CTRATrajectory;

} // namespace traji

/* helper functions for boost.geometry */
namespace boost { namespace geometry { namespace model
{
    inline bool operator==(const traji::Point& lhs, const traji::Point& rhs)
    {
        return lhs.get<0>() == rhs.get<0>() && lhs.get<1>() == rhs.get<1>();
    }
    inline bool operator!=(const traji::Point& lhs, const traji::Point& rhs)
    {
        return !(lhs == rhs);
    }
}}} // namespace boost.geometry.model

/* main implementations */
namespace traji
{

enum class SegmentType
{
    /// @brief A straight segment
    /// params: none
    Line,
    /// @brief A segment on a circle
    /// params: x, y of the arc center 
    Arc, 
    /// @brief A segment on a conic section
    /// params: TBD
    Conic,
    /// @brief A segment defined by a quadratic bezier curve
    /// params: x, y of the control point
    /// Note that a quadratic bezier curve is also a parabolic curve
    QuadraticBezier,
    /// @brief A segment defined by a cubic bezier curve
    /// params: x1, y1, x2, y2 of the two control points
    CubicBezier,

    // XXX: with all segment types implemented above, we should be able to export the path to a SVG

    /// @brief A path segment defined by Polynomial. This can be used to create a spline.
    /// params: polynomial coeffs (high to low)
    /// The latent parameter range need to be normalized to [0, 1].
    Polynomial, 
};

// ==================================== Non-parametric paths ====================================

/// Representation of a position in the path (and possibly extended line)
/// When the position is on the path, 0 <= segment < #segments, 0 <= fraction < 1
/// When the position is on the extended last segment, segment = #segments - 1, fraction >= 1
/// When the position is on the extended first segment, segment = 0, fraction < 0
// TODO: implement comparison operators
struct PathPosition
{
    std::size_t segment; // segment index in the line string
    TRel fraction; // fraction of the point in the segment

    PathPosition() : segment(0), fraction(0.0) {}
    PathPosition(std::size_t segment_, TRel fraction_)
        : segment(segment_), fraction(fraction_) {}

    /// Convert the position to the distance to the beginning
    TAbs to_s(const Path &path) const;
    TAbs to_s(const HeteroPath &path) const;

    /// Convert the distance to the beginning to s
    static PathPosition from_s(const Path &path, TAbs s);
    static PathPosition from_s(const HeteroPath &path, TAbs s);
    static std::vector<PathPosition> from_s(const Path &path, const std::vector<TAbs> &s_list);

    /// Convert the position to the distance to the timestamp
    TAbs to_t(const Trajectory &traj) const;
    static PathPosition from_t(const Trajectory &traj, TAbs t);
    static std::vector<PathPosition> from_t(const Trajectory &path, const std::vector<TAbs> &t_list);

    /// Move the position forward along the path. Note that this function is only efficient
    /// when the distance is not too large. Otherwise please use to_s() and from_s()
    PathPosition forward(const Path &path, TRel s) const;
    PathPosition backward(const Path &path, TRel s) const;
private:
    TRel to_rel_s(const Path &path) const;
    TRel to_rel_t(const Trajectory &traj) const;
    static PathPosition from_rel_s(const Path &path, TRel s);
    static PathPosition from_rel_t(const Trajectory &traj, TRel t);
};

/// (immutable) non-parametric linestring
class Path
{
public:
    /// @brief Distance of the first point along a certain path, default to 0. All s values will be offset by this amount when querying.
    TAbs s0 = 0;

protected:
    // The geometry of the line string
    LineString _line;

    // The distance of each vertices )starting from the second) to the first vertex along the line.
    // It has a length of n, where n is the number of vetices. The first value must be 0 (if exists).
    std::vector<TRel> _distance;

    /// Recalculate the _distance values
    void update_distance();

    /// Calculate feasible smooth radius given the input as max limit
    std::vector<TRel> calc_feasible_radius(TRel smooth_radius) const;

private:
    /// respacing with smooth_radius = 0
    Path respacing0(TRel resolution) const;
    TRel solve_curvature(size_t segment_idx) const;

public:
    friend class PathPosition;
    friend class QuinticPolyTrajectory;
    friend class CTRATrajectory;
    friend class HeteroPath;

    inline Path() {}
    template<typename Iterator>
    inline Path(Iterator begin, Iterator end, TRel s0_ = 0) : _line(begin, end), s0(s0_) { update_distance(); }
    inline Path(std::initializer_list<Point> l) : Path(l.begin(), l.end()) {}
    inline Path(const std::vector<Point> &l, TRel s0_ = 0): Path(l.begin(), l.end(), s0_) {}
    inline Path(std::vector<Point> &&l, TRel s0_ = 0): Path(l.begin(), l.end(), s0_) {}

    /// The size of a path is the number of segments
    inline std::size_t size() const { return _line.empty() ? 0 : _line.size() - 1; }
    inline TRel length() const { return _distance.back() - _distance.front(); }
    inline bool empty() const { return _line.empty(); }
    inline std::vector<TRel> segment_lengths() const
    {
        std::vector<TRel> lengths(size());
        for (int i = 1; i < _line.size(); i++)
            lengths[i-1] = _distance[i] - _distance[i-1];
        return lengths;
    }

    /// Access the underlying boost linestring
    inline LineString& data() { return _line; }
    inline const LineString& data() const { return _line; }

    /// Access vertices of the path
    inline std::vector<Point>& vertices() { return _line; }
    inline const LineString& vertices() const { return _line; }

    /// Get the point indicated by the distance from the beginning
    inline Point point_from(TAbs s) const
    {
        return point_at(PathPosition::from_s(*this, s));
    }
    Point point_at(const PathPosition &pos) const;

    /// Get the tagent at the point indicated by the distance from the beginning
    inline TRel tangent_from(TAbs s) const
    {
        return tangent_at(PathPosition::from_s(*this, s));
    }
    TRel tangent_at(const PathPosition &pos) const;

    /// Get the curvature at given position. The curvature is precise on the vertices, and
    /// interpolated on segments. The curvature is positive if the center is to the left
    /// of the path.
    inline TRel curvature_from(TAbs s) const
    {
        return curvature_at(PathPosition::from_s(*this, s));
    }
    TRel curvature_at(const PathPosition &pos) const;

    /// Get the value interpolated by the distance from the beginning
    inline TRel interpolate_from(const std::vector<TRel> &values, TAbs s) const
    {
        return interpolate_at(values, PathPosition::from_s(*this, s));
    }
    TRel interpolate_at(const std::vector<TRel> &values, const PathPosition &pos) const;

    /// Get the signed distance and foot point from a point to the path. The distance is positive
    /// if the point is at left hand side of the direction of path
    /// @return (distance to the path, projection point position)
    std::pair<TRel, PathPosition> project(const Point &point) const;

    // TODO: add a project_ext, which projects the point onto the extended line at the last segment.
    //       it should also allow specifying the expected path position. If the path contains loop,
    //       then we might just want a local minimum based on a starting s value

    /// Return the path in a form with equally space points (interpolating over s)
    /// The starting and end point of the original line string is guaranteed to be included and not smoothed
    /// @param smooth_radius The distance from original vertex to the new smoothed vertex.
    Path respacing(TRel resolution, TRel smooth_radius = 0) const;

    /// Return the path with each segment's length less than the given resolution
    /// This method will retain the original corner points, unlike respacing function
    Path densify(TRel resolution) const;

    /// Return the path with rounded corners. This method doesn't change the density of
    /// the points on the line segments
    HeteroPath smooth(TRel smooth_radius, SegmentType segtype = SegmentType::Arc) const;

    /// Return a path represented by a list of distance to start. 
    /// The best performance is achieved when s_list is sorted ascendingly.
    Path resample_from(const std::vector<TAbs> &s_list) const;

    /// Extract a part of the path (between s_start and s_end).
    /// The result path will contains all the way points between s_start and s_end,
    /// and the two end points (with extrapolating)
    Path extract_from(TAbs s_start, TAbs s_end) const;
};

/// non-parametric linestring with time, assuming constant acceleration
/// the values will still be interpolated based on position
class Trajectory : public Path
{
public:
    /// @brief Timestamp of the first point, default to 0. All t values will be offset by this amount when querying.
    TAbs t0;

protected:
    // The time elapsed from t0 at each vertices, ie timestamps but subtracted by t0.
    // It has a length of n, where n is the number of vetices. The first value must be 0.
    std::vector<TRel> _durations;
    
    // Calculate the duration values
    void update_duration(std::vector<TAbs> &&timestamps);

private:
    Vector2 solve_velocity(size_t segment_idx) const;
    Vector2 solve_acceleration(size_t point_idx) const;

public:
    friend class PathPosition;
    friend class QuinticPolyTrajectory;
    friend class CTRATrajectory;

    inline Trajectory() {}
    inline Trajectory(const Path& path, const std::vector<TAbs> &timestamps) : Path(path), _durations() {
        update_duration(std::vector<TAbs>(timestamps));
    }
    inline Trajectory(Path&& path, std::vector<TAbs> &&timestamps) : Path(path), _durations() {
        update_duration(std::vector<TAbs>(timestamps));
    }
    template<typename Iterator, typename TIterator>
    inline Trajectory(Iterator begin, Iterator end, TIterator t_begin, TIterator t_end, TAbs s0 = 0)
        : Path(begin, end, s0), _durations() {
        update_duration(std::vector<TAbs>(t_begin, t_end));
    }

    /// @brief Create the trajectory based on a path, assuming that the speed is constant along the trajectory.
    /// @param path The reference path
    /// @param t0 The initial timestamp
    /// @param speed The speed in the trajectory
    inline Trajectory(const Path& path, TAbs t0, TRel speed) : Path(path) {
        _durations.reserve(_distance.size());
        float d0 = _distance.front();
        for (auto &d : _distance) {
            _durations.push_back(t0 + (d - d0) / speed);
        }
    }

    std::vector<TAbs> timestamps() const;
    inline const TRel duration() const { return _durations.back() - _durations.front(); }

    /// Get the point indicated by the time
    using Path::point_at;
    inline Point point_at(TRel t) const 
    {
        return point_at(PathPosition::from_t(*this, t));
    }
    using Path::tangent_at;
    inline TRel tangent_at(TRel t) const
    {
        return tangent_at(PathPosition::from_t(*this, t));
    }
    Vector2 velocity_at(const PathPosition &pos, bool interpolate = false) const;
    inline Vector2 velocity_from(TRel s, bool interpolate = false) const
    {
        return velocity_at(PathPosition::from_s(*this, s), interpolate);
    }
    inline Vector2 velocity_at(TRel t, bool interpolate = false) const
    {
        return velocity_at(PathPosition::from_t(*this, t), interpolate);
    }
    Vector2 acceleration_at(const PathPosition &pos, bool interpolate = false) const;
    inline Vector2 acceleration_from(TRel s, bool interpolate = false) const
    {
        return acceleration_at(PathPosition::from_s(*this, s), interpolate);
    }
    inline Vector2 acceleration_at(TRel t, bool interpolate = false) const
    {
        return acceleration_at(PathPosition::from_t(*this, t), interpolate);
    }

    /// Return the trajectory represented by points with equal time intervals
    Trajectory reperiodize(TRel interval) const;

    /// Return a path represented by a list of timestamp. 
    /// The best performance is achieved when t_list is sorted ascendingly.
    Trajectory resample_at(const std::vector<TAbs> &t_list) const;
};


// ==================================== Parametric paths ====================================

/// This class represents an optimal trajectory determined by planning horizon T, starting states (s, v, a)
/// and end states at T. It's represented by a quintic polynomial. Optionally the final position on the x
/// direction can by relaxed by specify `relax_sx` on construction, then the optiimal trajectory is a quartic polynomial.
class QuinticPolyTrajectory
{
protected:
    Vector6 _x_coeffs, _y_coeffs; // high-order to low-order
    TRel _T; // the trajectory starts at t=0 and ends at t=T

public:
    QuinticPolyTrajectory(TRel T, const Vector6 &x_coeffs, const Vector6 &y_coeffs)
        : _x_coeffs(x_coeffs), _y_coeffs(y_coeffs), _T(T) {}

    /// Solve optimal trajectory based on x and y states (including 1st and 2nd state derivatives)
    QuinticPolyTrajectory(TRel T,
                          const Vector3 &x0, const Vector3 &xT,
                          const Vector3 &y0, const Vector3 &yT, bool relax_sx);

    inline const Vector6& x_coeffs() const { return _x_coeffs; }
    inline const Vector6& y_coeffs() const { return _y_coeffs; }
    inline TRel T() const { return _T; }

    Point point_at(TRel t) const;
    TRel tangent_at(TRel t) const;
    Vector2 velocity_at(TRel t) const;
    Vector2 acceleration_at(TRel t) const;

    // Trajectory rasterize(TRel resolution) const;
    Trajectory periodize(TRel interval) const;
};


/// Trajectory of a Constant Turn-Rate and (longitudinal) Acceleration model.
class CTRATrajectory
{
protected:
    Vector6 _init_state; // [x, y, theta, v, a, omega]
    TRel _T;

public:
    CTRATrajectory(TRel T, const Vector6 &init_state) : _T(T), _init_state(init_state) {}
    CTRATrajectory(TRel T, Point p, TRel theta, TRel v, TRel a = 0, TRel omega = 0)
        : _T(T), _init_state() {
        _init_state << p.get<0>(), p.get<1>(), theta, v, a, omega;
    }

    inline const Vector6& initial_state() const { return _init_state; }
    inline TRel T() const { return _T; }

    Point point_at(TRel t) const;
    inline TRel tangent_at(TRel t) const
    {
        return /*theta*/ _init_state(2) + /*omega*/ _init_state(5) * t;
    }
    Vector2 velocity_at(TRel t) const;

    Trajectory periodize(TRel interval) const;
};

// ==================================== Hybrid paths ====================================

struct HeteroSegment
{
    /// Type of the segment
    SegmentType type;

    /// Params for the segments, see SegmentType for details
    std::vector<TRel> params;

    TRel length (const Point &start, const Point &end) const;
    Point point_at (const Point &start, const Point &end, TRel fraction) const;
    TRel tangent_at (const Point &start, const Point &end, TRel fraction) const;
};

class HeteroPath
{
protected:
    std::vector<Point> _points;
    std::vector<HeteroSegment> _segments;
    std::vector<TRel> _distance; // same as Path::_distance

    /// Update _distance values [Deprecated interface]
    void update_distance(TRel s0 = 0);

public:
    HeteroPath() {}
    inline HeteroPath(const Path& path) : _points(path._line), _distance(path._distance),
        _segments(path.size(), HeteroSegment { SegmentType::Line, {} }) {}
    inline HeteroPath(Path&& path) : _points(std::move(path._line)), _distance(std::move(path._distance)),
        _segments(path.size(), HeteroSegment { SegmentType::Line, {} }) {}
    inline HeteroPath(const std::vector<Point> &points, const std::vector<HeteroSegment> &segments, TRel s0 = 0) :
        _points(points), _segments(segments) { update_distance(); }

    friend class PathPosition;
    friend class Path;

    /// Access vertices of the path
    inline std::vector<Point>& vertices() { return _points; }
    inline const std::vector<Point>& vertices() const { return _points; }

    /// number of segments
    inline std::size_t size() const { return _segments.size(); }
    /// total length of the hetero path
    inline TRel length() const { return _distance.back(); }
    /// size of each segment
    inline std::vector<TRel> segment_lengths() const
    {
        std::vector<TRel> lengths(size());
        for (int i = 1; i < _points.size(); i++)
            lengths[i-1] = _distance[i] - _distance[i-1];
        return lengths;
    }

    Point point_from(TAbs s)
    {
        return point_at(PathPosition::from_s(*this, s));
    }
    Point point_at(const PathPosition &pos);
    TRel tangent_from(TAbs s)
    {
        return tangent_at(PathPosition::from_s(*this, s));
    }
    TRel tangent_at(const PathPosition &pos);

    /// Convert to Path by only rasterizing curves (excl. line segment)
    /// If resolution <= 0, then the result will be directly generated from vertices
    Path rasterize(TRel resolution = 0) const;

    /// Return the path with each segment's length less than the given resolution
    HeteroPath densify(TRel resolution) const;
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
        return Trajectory(from_cartesian(ref, (Path&)traj), traj.timestamps());
    }
    inline Trajectory to_cartesian(const Path &ref, const Trajectory &traj)
    {
        return Trajectory(to_cartesian(ref, (Path&)traj), traj.timestamps());
    }

    // Convert the trajectory to cartesian coordinate and adjust vertices position to
    // accomodate curvature change
    Trajectory to_cartesian_rescale(const Path &ref, const Trajectory &traj);
}

// ==================================== binary functions ====================================

inline TRel distance(const Point &lhs, const Point &rhs) { return boost::geometry::distance(lhs, rhs); }
inline TRel distance(const Path &lhs, const Point &rhs) { return std::abs(lhs.project(rhs).first); }
inline PathPosition arg_distance(const Path &lhs, const Point &rhs) { return lhs.project(rhs).second; }

/// Return the closest distance at the same time point between two trajectories
TRel tdistance(const Trajectory &lhs, const Trajectory &rhs);

std::vector<Point> intersection(const Path &lhs, const Path &rhs);
// TODO: make the intersection result sorted?
std::vector<std::pair<PathPosition, PathPosition>> arg_intersection(const Path &lhs, const Path &rhs);

// ==================================== standard functions ====================================

inline bool operator==(const PathPosition& lhs, const PathPosition& rhs)
{
    return lhs.segment == rhs.segment && lhs.fraction == rhs.fraction;
}
inline bool operator!=(const PathPosition& lhs, const PathPosition& rhs) { return !(lhs == rhs); }
inline bool operator< (const PathPosition& lhs, const PathPosition& rhs)
{
    return lhs.segment == rhs.segment && lhs.fraction == rhs.fraction;
}
inline bool operator> (const PathPosition& lhs, const PathPosition& rhs) { return rhs < lhs; }
inline bool operator<=(const PathPosition& lhs, const PathPosition& rhs) { return !(lhs > rhs); }
inline bool operator>=(const PathPosition& lhs, const PathPosition& rhs) { return !(lhs < rhs); }

bool operator==(const Path& lhs, const Path& rhs);
inline bool operator!=(const Path& lhs, const Path& rhs) { return !(lhs == rhs); }

} // namespace traji

namespace std
{

std::string to_string(const traji::Point &value);
std::string to_string(const traji::Path &value);
std::string to_string(const traji::PathPosition &value);
std::string to_string(const traji::Trajectory &value);
std::string to_string(const traji::QuinticPolyTrajectory &value);

}

#endif // TRAJI_HPP
