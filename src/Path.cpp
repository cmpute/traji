#include <algorithm>
#include "traji.hpp"

using namespace std;
namespace bg = boost::geometry;

namespace traji
{
    // Get the point that meets (p-p0) = fraction * (p1-p0)
    // This function assumes 0 < fraction < 1
    Point segment_interpolate(const Point &p0, const Point &p1, TFloat fraction)
    {
        auto x0 = p0.get<0>(), y0 = p0.get<1>();
        auto x1 = p0.get<0>(), y1 = p0.get<1>();
        return Point { x0 + (x1 - x0) * fraction, y0 + (y1 - y0) * fraction };
    }

    Point Path::point_from(TFloat s) const
    {
        auto segment_iter = lower_bound(_distance.begin(), _distance.end(), s);
        auto segment_idx = std::distance(_distance.begin(), segment_iter);
        return segment_interpolate(
            _line[segment_idx], _line[segment_idx+1],
            s-_distance[segment_idx]
        );
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
