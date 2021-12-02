#include <cmath>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <string>
#include <iostream>
#include "traji.hpp"

using namespace std;
namespace bg = boost::geometry;

namespace traji
{
    /// Some helper functions for line segments
    namespace segment
    {
        /// Get the point that meets (p-p0) = fraction * (p1-p0)
        inline Point interpolate(const Point &p0, const Point &p1, TFloat fraction)
        {
            auto x0 = p0.get<0>(), y0 = p0.get<1>();
            auto x1 = p1.get<0>(), y1 = p1.get<1>();
            return Point { x0 + (x1 - x0) * fraction, y0 + (y1 - y0) * fraction };
        }

        /// Calculate the signed distance from point to the segment **line**
        /// @return (distance to line, fraction of the foot point)
        /// The distance is positive if the point is at left hand side of the direction of line (p0 -> p1)
        pair<TFloat, TFloat> sdistance(const Point &p0, const Point &p1, const TFloat l, const Point &p)
        {
            auto x0 = p0.get<0>(), y0 = p0.get<1>();
            auto x1 = p1.get<0>(), y1 = p1.get<1>();
            auto x = p.get<0>(), y = p.get<1>();

            auto dx = x1 - x0, dy = y1 - y0;

            // if two points are the same one
            if (l == 0)
            {
                TFloat ds = hypot(x - x0, y - y0);
                return make_pair(ds, 0);
            }

            auto ds = (dx*y - dy*x + x0*y1 - x1*y0) / l;
            auto d0 = (x0*x0+x*dx-x0*x1 + y0*y0+y*dy-y0*y1) / l; // distance from foot point to p0
            
            return make_pair(ds, d0 / l);
        }

        inline TFloat tangent(const Point &p0, const Point &p1)
        {
            return atan2(p1.get<1>() - p0.get<1>(), p1.get<0>() - p0.get<0>());
        }

        /// Convert the sdistance to normal unsigned (squared) distance
        inline TFloat distance2(const pair<TFloat, TFloat> &sdist, const TFloat l)
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

    // Some helper function for arc segments
    namespace arc
    {
        struct ArcParams
        {
            Point center;
            TFloat radius;
            TFloat angle; // angle < 0 means center is on the right of the line
            TFloat start_angle;
        };

        /// Wrap the input angle to -pi~pi range
        inline TFloat warp_angle(TFloat angle)
        {
            return fmod(angle + pi, 2 * pi) - pi;
        }

        /// Solve the parameters of the smoothing arc between two adjacent segments
        ArcParams solve_smooth(
            const Point &p0, const Point &pivot, const Point &p1,
            TFloat l0, TFloat l1, TFloat smooth_radius)
        {
            auto tan1 = segment::tangent(p0, pivot);
            auto tan2 = segment::tangent(pivot, p1);

            auto turn_angle = warp_angle(tan2 - tan1);
            auto half = abs(turn_angle / 2);
            auto p2c_dist = smooth_radius / cos(half);

            TFloat radius = smooth_radius * tan(half);
            TFloat angle_start = tan1 + (turn_angle > 0 ? -pi2 : pi2);
            TFloat angle_mid = angle_start + (turn_angle > 0 ? half : -half);

            Point center(pivot.get<0>() - p2c_dist * cos(angle_mid),
                         pivot.get<1>() - p2c_dist * sin(angle_mid));

            TFloat angle = pi - abs(turn_angle);

            return ArcParams {center, radius, turn_angle > 0 ? angle : -angle, angle_start};
        }

        inline Point interpolate(const ArcParams &params, TFloat fraction)
        {
            auto angular_pos = params.start_angle + params.angle * fraction;

            return Point(
                params.center.get<0>() + cos(angular_pos) * params.radius,
                params.center.get<1>() + sin(angular_pos) * params.radius
            );
        }
    }

    TFloat PathPosition::to_s(const Path &path)
    {
        return path._distance[segment] + (path._distance[segment+1] - path._distance[segment]) * fraction;
    }

    PathPosition PathPosition::from_s(const Path &path, TFloat s)
    {
        auto segment_iter = lower_bound(path._distance.begin(), path._distance.end(), s);
        auto segment_idx = distance(path._distance.begin(), segment_iter);
        auto s0 = path._distance[segment_idx], s1 = path._distance[segment_idx+1];

        PathPosition result;
        result.segment = segment_idx;
        result.fraction = (s - s0) / (s1 - s0);
        return result;
    }

    void Path::update_distance()
    {
        if (_line.size() == 1) // Fix invalid line string (only 1 point)
            _line.clear();

        TFloat s = 0;
        _distance.resize(_line.size());
        _distance[0] = 0;
        for (int i = 1; i < _line.size(); i++)
        {
            s += distance(_line[i], _line[i-1]);
            _distance[i] = s;
        }
    }

    vector<PathPosition> PathPosition::from_s(const Path &path, const vector<TFloat> &s_list)
    {
        vector<PathPosition> result;
        result.reserve(s_list.size());

        TFloat cur_s = 0;
        size_t cur_idx = 1; // index of first point that has s larger than cur_s
        for (auto s : s_list)
            if (s < cur_s)
                result.push_back(PathPosition::from_s(path, s));
            else
            {
                while (s > path._distance[cur_idx] && cur_idx < (path.size() - 1))
                    cur_idx++;

                PathPosition pos;
                pos.segment = cur_idx - 1;
                auto s0 = path._distance[cur_idx-1], s1 = path._distance[cur_idx];
                pos.fraction = (s - s0) / (s1 - s0);
                result.push_back(pos);

                cur_s = s;
            }

        return result;
    }

    Point Path::point_at(const PathPosition &pos) const
    {
        return segment::interpolate(
            _line[pos.segment], _line[pos.segment+1],
            pos.fraction
        );
    }

    TFloat Path::tangent_at(const PathPosition &pos) const
    {
        return segment::tangent(
            _line[pos.segment], _line[pos.segment+1]
        );
    }

    pair<TFloat, PathPosition> Path::project(const Point &point) const
    {
        vector<pair<TFloat, TFloat>> dists(_line.size() - 1);
        vector<TFloat> comp_dists(_line.size() - 1); // actual distance for comparison
        for (int i = 1; i < _line.size(); i++)
        {
            // calculate signed distance
            int iprev = i - 1;
            auto l = _distance[i] - _distance[iprev];
            dists[iprev] = segment::sdistance(_line[iprev], _line[i], l, point);

            // calculate squared distance for comparison
            comp_dists[iprev] = segment::distance2(dists[iprev], l);
        }

        auto min_idx = distance(comp_dists.begin(),
            min_element(comp_dists.begin(), comp_dists.end()));

        pair<TFloat, PathPosition> result;
        result.first = dists[min_idx].first;
        result.second.segment = min_idx;
        result.second.fraction = dists[min_idx].second;
        return result;
    }

    Path Path::densify(TFloat resolution) const
    {
        Path result;
        result._line.push_back(_line.front());
        for (size_t i = 1; i < _line.size(); i++)
        {
            int mul = ceil((_distance[i] - _distance[i-1]) / resolution);
            if (mul <= 1)
                result._line.push_back(_line[i]);
            else
            {
                for (size_t j = 1; j <= mul; j++)
                    result._line.push_back(segment::interpolate(_line[i-1], _line[i], (TFloat) j / mul));
            }
        }

        result.update_distance();
        return result;
    }

    vector<TFloat> Path::calc_feasible_radius(TFloat smooth_radius) const
    {
        assert(_line.size() > 2);
        vector<TFloat> radius(_line.size() - 2, smooth_radius);
        vector<TFloat> seglen = segment_lengths();

        // sort segment lengths
        vector<size_t> indices(seglen.size());
        iota(indices.begin(), indices.end(), 0);
        sort(indices.begin(), indices.end(),
            [&seglen](int left, int right) -> bool {
                // sort indices according to corresponding array element
                return seglen[left] < seglen[right];
            });

        // update maximum radius
        for (auto i : indices)
        {
            if (seglen[i] > 2 * smooth_radius)
                break;

            if (i == 0)
                radius[0] = min(radius[0], seglen[0]);
            else if (i == radius.size()) // n-2
                radius[i-1] = min(radius[i-1], seglen[i]);
            else
            {
                radius[i-1] = min(radius[i-1], seglen[i] / 2);
                radius[i] = min(radius[i], seglen[i] / 2);
            }
        }
        return radius;
    }

    Path Path::respacing0(TFloat resolution) const
    {
        assert (_line.size() >= 2);

        Path result;
        result._line.push_back(_line.front());

        auto residual_s = resolution;
        auto segment_length = _distance[1];
        size_t segment_idx = 0;
        
        while(residual_s < segment_length)
        {
            result._line.push_back(segment::interpolate(
                _line[segment_idx], _line[segment_idx+1],
                residual_s / segment_length
            ));

            residual_s += resolution;
            if (residual_s >= segment_length)
            {
                residual_s -= segment_length;

                segment_idx++;
                if (segment_idx >= _line.size() - 1)
                    break;
                segment_length = _distance[segment_idx+1] - _distance[segment_idx];
            }
        }

        if (residual_s > segment_length)
            result._line.push_back(_line.back());

        result.update_distance();
        return result;
    }

    Path Path::respacing(TFloat resolution, TFloat smooth_radius) const
    {
        // shortcuts
        if (_line.size() == 0)
            return *this;
        if (smooth_radius == 0 || _line.size() == 2)
            return respacing0(resolution);

        Path result;
        result._line.push_back(_line.front());
        vector<TFloat> seglen = segment_lengths();

        auto residual_s = resolution;
        vector<TFloat> feas_radius = calc_feasible_radius(smooth_radius);
        auto segment_length = seglen[0] - feas_radius[0]; // current effective segment length
        TFloat segment_soffset = 0; // start offset
        size_t segment_idx = 0;

        bool on_arc = false;
        arc::ArcParams arc_params;

        while (true)
        {
            if (on_arc)
            {
                if (residual_s >= segment_length)
                {
                    // jump to the next segment
                    residual_s -= segment_length;
                    on_arc = false;

                    segment_idx++;
                    if (segment_idx == _line.size() - 2)
                        segment_length = seglen[segment_idx] - feas_radius[segment_idx-1];
                    else
                        segment_length = seglen[segment_idx] - feas_radius[segment_idx-1] - feas_radius[segment_idx];
                    segment_soffset = feas_radius[segment_idx-1];
                    continue;
                }
                else
                {
                    // interpolate on arc
                    result._line.push_back(arc::interpolate(
                        arc_params,
                        residual_s / segment_length
                    ));
                    residual_s += resolution;
                }
            }
            else // on segment
            {
                if (residual_s >= segment_length)
                {
                    // jump to the next arc
                    if (segment_idx >= _line.size() - 2)
                        break;

                    residual_s -= segment_length;

                    on_arc = true;
                    arc_params = arc::solve_smooth(
                        _line[segment_idx], _line[segment_idx+1], _line[segment_idx+2],
                        seglen[segment_idx], seglen[segment_idx+1],
                        feas_radius[segment_idx]);
                    segment_length = abs(arc_params.radius * arc_params.angle);
                }
                else
                {
                    // interpolate on line segment
                    result._line.push_back(segment::interpolate(
                        _line[segment_idx], _line[segment_idx+1],
                        (residual_s + segment_soffset) / seglen[segment_idx]
                    ));
                    residual_s += resolution;
                }
            }
        }

        if (residual_s > segment_length)
            result._line.push_back(_line.back());

        result.update_distance();
        return result;
    }

    HeteroPath Path::smooth(TFloat smooth_radius) const
    {
        HeteroPath result;
        if (_line.size() == 0)
            return result;
        if (_line.size() == 2)
        {
            result._segments.emplace_back(HeteroSegment {
                SegmentType::Line, _line[0], _line[1], std::vector<TFloat>()
            });
            result._distance = {_distance[1]};
            return result;
        }

        auto feas_radius = calc_feasible_radius(smooth_radius);
        auto seglen = segment_lengths();
        result._segments.reserve(seglen.size());
        result._distance.reserve(seglen.size());

        // add first segment
        auto last_point = segment::interpolate(
            _line[0], _line[1],
            1 - feas_radius[0] / seglen[0]
        );
        result._segments.emplace_back(HeteroSegment {
            SegmentType::Line, _line.front(), last_point, std::vector<TFloat>()
        });
        result._distance.push_back(seglen.front() - feas_radius.front());

        for (size_t i = 0; i < seglen.size() - 1; i++)
        {
            // add i-th segment
            if (i > 0)
            {
                auto new_point = segment::interpolate(
                    _line[i], _line[i+1],
                    1 - feas_radius[i] / seglen[i]
                );
                result._segments.emplace_back(HeteroSegment {
                    SegmentType::Line, last_point, new_point,
                    std::vector<TFloat>()
                });
                result._distance.push_back(result._distance.back() +
                    seglen[i] - feas_radius[i-1] - feas_radius[i]);
                last_point = new_point;
            }

            // add i-th curve
            auto arc_params = arc::solve_smooth(
                _line[i], _line[i+1], _line[i+2],
                seglen[i], seglen[i+1],
                feas_radius[i]);
            auto arc_len = abs(arc_params.angle * arc_params.radius);

            auto new_point = segment::interpolate(
                _line[i+1], _line[i+2],
                feas_radius[i] / seglen[i+1]
            );
            result._segments.emplace_back(HeteroSegment {
                SegmentType::Arc, last_point, new_point,
                std::vector<TFloat> {arc_params.center.get<0>(), arc_params.center.get<1>()}
            });
            result._distance.push_back(result._distance.back() + arc_len);
            last_point = new_point;
        }

        // add last segment
        result._segments.emplace_back(HeteroSegment {
            SegmentType::Line, last_point, _line.back(),
            std::vector<TFloat>()
        });
        result._distance.push_back(seglen.back() - feas_radius.back());

        return result;
    }

    Path Path::resample(const vector<TFloat> &s_list) const
    {
        Path result;
        result._line.reserve(s_list.size());
        vector<PathPosition> pos_list = PathPosition::from_s(*this, s_list);
        for (auto pos : pos_list)
            result._line.push_back(point_at(pos));
    }

    vector<Point> intersection(const Path &lhs, const Path &rhs)
    {
        vector<Point> plist;
        bg::intersection(lhs.data(), rhs.data(), plist);
        return plist;
    }

    vector<pair<PathPosition, PathPosition>> arg_intersection(const Path &lhs, const Path &rhs)
    {
        // TODO: implement this arg intersection algorithm using sweeping line or even
        //       Bentley-Ottmann algorithm to avoid calling "project"
        // Also need to benchmark the performance vs boost version

        vector<Point> plist = intersection(lhs, rhs);
        vector<pair<PathPosition, PathPosition>> result (plist.size());
        std::transform(plist.begin(), plist.end(), result.begin(),
            [&lhs, &rhs](const Point& p){ return make_pair(lhs.project(p).second, rhs.project(p).second); });
        return result;
    }

    bool operator==(const Path& lhs, const Path& rhs)
    {
        if(lhs.data().size() != rhs.data().size())
            return false;
        for(size_t i = 0; i < lhs.data().size(); i++)
            if (lhs[i] != rhs[i])
                return false;
        return true;
    }
}

namespace std
{
    string to_string(const traji::Point &value)
    {
        stringstream ss;
        ss << "(" << value.get<0>() << ", " << value.get<1>() << ")";
        return ss.str();
    }

    string to_string(const traji::Path &value)
    {
        if (value.size() == 0)
            return string("[]");

        stringstream ss;
        ss << '[' << to_string(value[0]);
        for (size_t i = 1; i < value.size(); i++)
            ss << ", " << to_string(value[i]);
        ss << ']';
        return ss.str();
    }

    string to_string(const traji::PathPosition &pos)
    {
        stringstream ss;
        ss << "(component " << pos.component << ", segment " << pos.segment << ", fraction " << pos.fraction << ')';
        return ss.str();
    }
}
