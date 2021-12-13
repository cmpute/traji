#include <cmath>
#include <algorithm>
#include <numeric>
#include <string>
#include <iostream>
#include "traji.hpp"
#include "traji_impl.hpp"

using namespace std;
using namespace std::placeholders;
namespace bg = boost::geometry;

namespace traji
{
    TFloat PathPosition::to_s(const Path &path)
    {
        return path._distance[segment] + (path._distance[segment+1] - path._distance[segment]) * fraction;
    }

    PathPosition PathPosition::from_s(const Path &path, TFloat s)
    {
        auto segment_iter = upper_bound(path._distance.begin(), path._distance.end(), s);
        auto segment_idx = distance(path._distance.begin(), segment_iter) - 1;
        segment_idx = min((long)path.size() - 2L, max(0L, segment_idx)); // clip the segment index

        auto s0 = path._distance[segment_idx], s1 = path._distance[segment_idx+1];

        return PathPosition(segment_idx, (s - s0) / (s1 - s0));
    }

    void Path::update_distance(TFloat s0)
    {
        if (_line.size() == 1) // Fix invalid line string (only 1 point)
        {
            _line.clear();
            return;
        }

        TFloat s = s0;
        _distance.resize(_line.size());
        _distance[0] = s;
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

                auto s0 = path._distance[cur_idx-1], s1 = path._distance[cur_idx];
                result.push_back(PathPosition(cur_idx - 1, (s - s0) / (s1 - s0)));

                cur_s = s;
            }

        return result;
    }

    Point Path::point_at(const PathPosition &pos) const
    {
        return line::interpolate(
            _line[pos.segment], _line[pos.segment+1],
            pos.fraction
        );
    }

    TFloat Path::tangent_at(const PathPosition &pos) const
    {
        return line::tangent(
            _line[pos.segment], _line[pos.segment+1]
        );
    }

    TFloat Path::interpolate_at(const std::vector<TFloat> &values, const PathPosition &pos) const
    {
        assert (values.size() == size());
        return values[pos.segment] * pos.fraction + values[pos.segment+1] * (1 - pos.fraction);
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
            dists[iprev] = line::sdistance(_line[iprev], _line[i], l, point);

            // calculate squared distance for comparison
            comp_dists[iprev] = line::distance2(dists[iprev], l);
        }

        auto min_idx = distance(comp_dists.begin(),
            min_element(comp_dists.begin(), comp_dists.end()));

        return make_pair(dists[min_idx].first,
            PathPosition(min_idx, dists[min_idx].second));
    }

    Path Path::densify(TFloat resolution) const
    {
        if (size() == 0)
            return *this;

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
                    result._line.push_back(line::interpolate(_line[i-1], _line[i], (TFloat) j / mul));
            }
        }

        result.update_distance(_distance[0]);
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
            result._line.push_back(line::interpolate(
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

        result.update_distance(_distance[0]);
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
                    result._line.push_back(line::interpolate(
                        _line[segment_idx], _line[segment_idx+1],
                        (residual_s + segment_soffset) / seglen[segment_idx]
                    ));
                    residual_s += resolution;
                }
            }
        }

        if (residual_s > segment_length)
            result._line.push_back(_line.back());

        result.update_distance(_distance[0]);
        return result;
    }

    HeteroPath Path::smooth(TFloat smooth_radius, SegmentType segtype) const
    {
        HeteroPath result;

        // shortcuts
        if (_line.size() == 0)
            return result;
        if (_line.size() == 2)
        {
            result._points = { _line[0], _line[1] };
            result._segments.emplace_back(HeteroSegment {
                SegmentType::Line, std::vector<TFloat>()
            });
            result._distance = { _distance[0], _distance[1] };
            return result;
        }

        auto feas_radius = calc_feasible_radius(smooth_radius);
        auto seglen = segment_lengths();
        result._points.reserve(_line.size() * 2); // TODO: better estimate with correct Path::size()
        result._segments.reserve(seglen.size());
        result._distance.reserve(seglen.size());

        // add first segment
        result._points.push_back(_line[0]);
        result._points.emplace_back(line::interpolate(
            _line[0], _line[1],
            1 - feas_radius[0] / seglen[0]
        ));
        result._segments.emplace_back(HeteroSegment {
            SegmentType::Line, std::vector<TFloat>()
        });
        result._distance.emplace_back(seglen.front() - feas_radius.front());

        // loop over the intermediate segments
        for (size_t i = 0; i < seglen.size() - 1; i++)
        {
            // add i-th segment
            if (i > 0)
            {
                result._points.emplace_back(line::interpolate(
                    _line[i], _line[i+1],
                    1 - feas_radius[i] / seglen[i]
                ));
                result._segments.emplace_back(HeteroSegment { SegmentType::Line, {} });
                result._distance.emplace_back(result._distance.back() +
                    seglen[i] - feas_radius[i-1] - feas_radius[i]);
            }

            // add i-th curve
            result._points.emplace_back(line::interpolate(
                _line[i+1], _line[i+2],
                feas_radius[i] / seglen[i+1]
            ));
            switch (segtype)
            {
                case SegmentType::Line:
                {
                    result._segments.emplace_back(HeteroSegment { SegmentType::Line, {} });
                    result._distance.emplace_back(distance(
                        result._points.rbegin()[1], result._points.rbegin()[0]));
                    break;
                }
                case SegmentType::Arc:
                {
                    auto arc_params = arc::solve_smooth(
                        _line[i], _line[i+1], _line[i+2],
                        seglen[i], seglen[i+1],
                        feas_radius[i]);
                    auto arc_len = abs(arc_params.angle * arc_params.radius);

                    result._segments.emplace_back(HeteroSegment {
                        SegmentType::Arc,
                        { arc_params.center.get<0>(), arc_params.center.get<1>() }
                    });
                    result._distance.emplace_back(result._distance.back() + arc_len);
                    break;
                }
                default:
                    throw std::runtime_error("Unsupported segment type!");
            }
        }

        // add last segment
        result._points.push_back(_line.back());
        result._segments.emplace_back(HeteroSegment {
            SegmentType::Line, std::vector<TFloat>()
        });
        result._distance.emplace_back(seglen.back() - feas_radius.back());

        return result;
    }

    Path Path::resample_from(const vector<TFloat> &s_list) const
    {
        vector<PathPosition> pos_list = PathPosition::from_s(*this, s_list);
        vector<Point> plist; plist.reserve(s_list.size());
        transform(pos_list.begin(), pos_list.end(),
            std::back_inserter(plist), std::bind(&Path::point_at, this, _1));
        return Path(move(plist));
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
        vector<pair<PathPosition, PathPosition>> result; result.reserve(plist.size());
        transform(plist.begin(), plist.end(), back_inserter(result),
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
