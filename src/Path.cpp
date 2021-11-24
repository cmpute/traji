#include <cmath>
#include <algorithm>
#include "traji.hpp"

using namespace std;
namespace bg = boost::geometry;

namespace traji
{
    Point Segment::interpolate(const Point &p0, const Point &p1, TFloat fraction)
    {
        auto x0 = p0.get<0>(), y0 = p0.get<1>();
        auto x1 = p0.get<0>(), y1 = p0.get<1>();
        return Point { x0 + (x1 - x0) * fraction, y0 + (y1 - y0) * fraction };
    }

    std::pair<TFloat, TFloat> Segment::sdistance(
        const Point &p0, const Point &p1, const TFloat l, const Point &p)
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

    Point Path::point_from(TFloat s) const
    {
        auto segment_iter = lower_bound(_distance.begin(), _distance.end(), s);
        auto segment_idx = std::distance(_distance.begin(), segment_iter);
        auto d0 = _distance[segment_idx];
        return Segment::interpolate(
            _line[segment_idx], _line[segment_idx+1],
            (s - d0) / (_distance[segment_idx+1] - d0)
        );
    }

    Point Path::point_at(const PathPosition &pos) const
    {
        return Segment::interpolate(
            _line[pos.segment], _line[pos.segment+1],
            pos.fraction
        );
    }

    TFloat Path::tangent_from(TFloat s) const
    {
        // TODO: add round_radius param
        // TODO: use bg::azimuth
    }

    PathPosition Path::project(const Point &point) const
    {
        // TODO
    }

    std::pair<TFloat, PathPosition> Path::sdistance(const Point &point) const
    {
        std::vector<std::pair<TFloat, TFloat>> dists(_line.size() - 1);
        std::vector<TFloat> comp_dists(_line.size() - 1); // actual distance for comparison
        for (int i = 1; i < dists.size(); i++)
        {
            // calculate signed distance
            int iprev = i - 1;
            auto l = _distance[i] - _distance[iprev];
            dists[iprev] = Segment::sdistance(_line[iprev], _line[i], l, point);

            // calculate distance
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
            [&lhs, &rhs](const Point& p){ return make_pair(lhs.project(p), rhs.project(p)); });
        return result;
    }
}
