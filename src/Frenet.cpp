#include <algorithm>
#include <cmath>
#include "traji.hpp"

using namespace std;

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
    vector<Point> frenet_points(path.size());
    transform(path.data().begin(), path.data().end(), frenet_points.begin(),
        [&ref](auto p) { return from_cartesian(ref, p); });
    return Path(frenet_points.begin(), frenet_points.end());
}

Path to_cartesian(const Path &ref, const Path &path)
{
    vector<Point> cartesian_points(path.size());
    vector<TFloat> s_list(path.size()); // XXX: avoid copy if the underlying structure is changed to Eigen
    transform(path.data().begin(), path.data().end(), s_list.begin(),
              [](const Point &p) { return p.get<0>(); });
    vector<PathPosition> pos_list = PathPosition::from_s(ref, s_list);
    for (size_t i = 0; i < path.size(); i++)
        cartesian_points[i] = to_cartesian_hint(ref, path.data()[i], pos_list[i]);
    return Path(cartesian_points.begin(), cartesian_points.end());
}

} // namespace frenet    
} // namespace traji
