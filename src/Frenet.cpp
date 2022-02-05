#include <algorithm>
#include <numeric>
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

    auto pos_list = PathPosition::from_s(ref, s_list);
    vector<Point> cartesian_points; cartesian_points.reserve(path.size() + 1);
    for (size_t i = 0; i <= path.size(); i++)
        cartesian_points.push_back(to_cartesian_hint(ref, path.vertices()[i], pos_list[i]));
    return Path(move(cartesian_points));
}

Trajectory to_cartesian_fixing_position(const Path &ref, const Trajectory &traj, TFloat scale)
{
    // Calculate the original path without resampling
    vector<TFloat> s_list; s_list.reserve(traj.size() + 1);
    transform(traj.vertices().begin(), traj.vertices().end(), back_inserter(s_list),
              [](const Point &p) { return p.get<0>(); });

    auto pos_list = PathPosition::from_s(ref, s_list);
    vector<Point> cartesian_points; cartesian_points.reserve(traj.size() + 1);
    for (size_t i = 0; i <= traj.size(); i++)
        cartesian_points.push_back(to_cartesian_hint(ref, traj.vertices()[i], pos_list[i]));
    auto original_path = Path(move(cartesian_points));

    // Adjust segment lengths by curv_arr
    ArrayX n_list(traj.size() + 1);
    for (size_t i = 0; i <= traj.size(); i++)
        n_list(i) = traj.vertices()[i].get<1>();

    auto seglen = traj.segment_lengths();
    auto segcount = seglen.size();
    ArrayX curv_arr(segcount);
    for (size_t i = 0; i <= traj.size(); i++)
        curv_arr[i] = ref.curvature_at(pos_list[i]);

    auto seglen_arr = ArrayX::Map(&seglen[0], segcount);
    auto scale_factor = (1 - n_list * curv_arr) * scale;
    auto scale_normed = scale_factor / scale_factor.sum();
    vector<TFloat> scaled_seglen; scaled_seglen.reserve(segcount+1);
    scaled_seglen.push_back(0); // this is s0 for the new path
    for (size_t i = 0; i < segcount; i++)
        scaled_seglen.push_back(seglen_arr[i] / scale_normed[i]);
    vector<TFloat> new_s_list(segcount+1);
    partial_sum(scaled_seglen.begin(), scaled_seglen.end(), new_s_list.begin(), plus<TFloat>());
    
    // Calculate new points
    return Trajectory(original_path.resample_from(new_s_list), traj.timestamps());
}

} // namespace frenet    
} // namespace traji
