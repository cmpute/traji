#include <algorithm>
#include <cmath>
#include "traji.hpp"

using namespace std;
using namespace std::placeholders;

namespace traji { namespace frenet
{

Point from_cartesian(const Path &ref, const Point &point)
{
    auto result = ref.project(point);
    auto distance = result.first;
    auto proj = result.second;
    return Point(proj.to_s(ref), distance);
}

inline Point to_cartesian_hint(const Path &ref, const Point &point, const PathPosition &s_pos)
{
    auto tan = ref.tangent_at(s_pos);
    auto base_point = ref.point_at(s_pos);
    auto normal = tan + pi2;
    auto l = point.get<1>();
    return Point(base_point.get<0>() + l * cos(normal),
                 base_point.get<1>() + l * sin(normal));
}

Point to_cartesian(const Path &ref, const Point &point)
{
    return to_cartesian_hint(ref, point, PathPosition::from_s(ref, point.get<0>()));
}

Path from_cartesian(const Path &ref, const Path &path)
{
    vector<Point> frenet_points; frenet_points.reserve(path.size() + 1);
    transform(path.vertices().begin(), path.vertices().end(),
        back_inserter(frenet_points),
        bind<Point(const Path&, const Point&)>(from_cartesian, ref, _1));
    return Path(move(frenet_points));
}

Path to_cartesian(const Path &ref, const Path &path)
{
    vector<TFloat> s_list; s_list.reserve(path.size() + 1);
    transform(path.vertices().begin(), path.vertices().end(), back_inserter(s_list),
              [](const Point &p) { return p.get<0>(); });

    vector<PathPosition> pos_list = PathPosition::from_s(ref, s_list);
    vector<Point> cartesian_points; cartesian_points.reserve(path.size() + 1);
    for (size_t i = 0; i <= path.size(); i++)
        cartesian_points.push_back(to_cartesian_hint(ref, path.vertices()[i], pos_list[i]));
    return Path(move(cartesian_points));
}

} // namespace frenet    
} // namespace traji
