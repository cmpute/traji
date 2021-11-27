#include <cmath>
#include <algorithm>
#include <string>
#include "traji.hpp"

using namespace std;
namespace bg = boost::geometry;

namespace traji
{
    /// Some helper functions for line segments
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
        static std::pair<TFloat, TFloat> sdistance(const Point &p0, const Point &p1, const TFloat l, const Point &p)
        {
            auto x0 = p0.get<0>(), y0 = p0.get<0>();
            auto x1 = p1.get<0>(), y1 = p1.get<1>();
            auto x = p.get<0>(), y = p.get<1>();

            auto dx = x1 - x0, dy = y1 - y0;

            // if two points are the same one
            if (l == 0)
            {
                TFloat ds = hypot(x - x0, y - y0);
                return std::make_pair(ds, 0);
            }

            auto ds = (dx*y - dy*x + x0*y1 - x1*y0) / l;
            auto d0 = (x0*x0+x*dx-x0*x1 + y0*y0+y*dy-y0*y1) / l; // distance from foot point to p0
            
            return make_pair(ds, d0 / l);
        }

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

    TFloat PathPosition::to_s(const Path &path)
    {
        return path._distance[segment] + (path._distance[segment+1] - path._distance[segment]) * fraction;
    }

    PathPosition PathPosition::from_s(const Path &path, TFloat s)
    {
        auto segment_iter = lower_bound(path._distance.begin(), path._distance.end(), s);
        auto segment_idx = std::distance(path._distance.begin(), segment_iter);
        auto s0 = path._distance[segment_idx];

        PathPosition result;
        result.segment = segment_idx;
        result.fraction = (s - s0) / (path._distance[segment_idx+1] - s0);
        return result;
    }

    void Path::update_distance()
    {
        TFloat s = 0;
        _distance.resize(_line.size());
        _distance[0] = 0;
        for (int i = 1; i < _line.size(); i++)
        {
            s += distance(_line[i], _line[i-1]);
            _distance[i] = s;
        }
    }

    Point Path::point_at(const PathPosition &pos) const
    {
        return Segment::interpolate(
            _line[pos.segment], _line[pos.segment+1],
            pos.fraction
        );
    }

    TFloat Path::tangent_at(const PathPosition &pos) const
    {
        auto p0 = _line[pos.segment+1], p1 = _line[pos.segment];
        auto segment_tan = std::atan2(p1.get<1>() - p0.get<1>(), p1.get<0>() - p0.get<0>());
        return segment_tan;
    }

    std::pair<TFloat, PathPosition> Path::project(const Point &point) const
    {
        std::vector<std::pair<TFloat, TFloat>> dists(_line.size() - 1);
        std::vector<TFloat> comp_dists(_line.size() - 1); // actual distance for comparison
        for (int i = 1; i < dists.size(); i++)
        {
            // calculate signed distance
            int iprev = i - 1;
            auto l = _distance[i] - _distance[iprev];
            dists[iprev] = Segment::sdistance(_line[iprev], _line[i], l, point);

            // calculate squared distance for comparison
            comp_dists[iprev] = Segment::distance2(dists[iprev], l);
        }

        auto min_idx = std::distance(comp_dists.begin(),
            std::min_element(comp_dists.begin(), comp_dists.end()));

        std::pair<TFloat, PathPosition> result;
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
            int mul = std::ceil((_distance[i] - _distance[i-1]) / resolution);
            if (mul <= 1)
                result._line.push_back(_line[i]);
            else
            {
                for (int i = 0; i < mul; i++)
                    result._line.push_back(Segment::interpolate(_line[i-1], _line[i], (TFloat) i / mul));
            }
        }

        result.update_distance();
        return result;
    }

    Path Path::respacing(TFloat resolution, float smooth_radius) const
    {
        if (smooth_radius != 0)
            throw "Not implemented";

        Path result;
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
        //       Bentley–Ottmann algorithm to avoid calling "project"
        // Also need to benchmark the performance vs boost version

        vector<Point> plist = intersection(lhs, rhs);
        vector<pair<PathPosition, PathPosition>> result (plist.size());
        std::transform(plist.begin(), plist.end(), result.begin(),
            [&lhs, &rhs](const Point& p){ return make_pair(lhs.project(p).second, rhs.project(p).second); });
        return result;
    }
}

namespace std
{
    std::string to_string(const traji::Point &value)
    {
        std::stringstream ss;
        ss << "(" << value.get<0>() << ", " << value.get<1>() << ")";
        return ss.str();
    }

    std::string to_string(const traji::Path &value)
    {
        if (value.size() == 0)
            return std::string("[]");

        std::stringstream ss;
        ss << '[' << to_string(value[0]);
        for (size_t i = 1; i < value.size(); i++)
            ss << ", " << to_string(value[i]);
        ss << ']';
        return ss.str();
    }

    std::string to_string(const traji::PathPosition &pos)
    {
        std::stringstream ss;
        ss << "(component " << pos.component << ", segment " << pos.segment << ", fraction" << pos.fraction << ')';
        return ss.str();
    }
}
