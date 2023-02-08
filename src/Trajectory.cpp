#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "traji.hpp"

using namespace std;
using namespace std::placeholders;

namespace traji
{
    TAbs PathPosition::to_t(const Trajectory &traj) const
    {
        return traj._durations[segment] + (traj._durations[segment+1] - traj._durations[segment]) * fraction;
    }
    PathPosition PathPosition::from_t(const Trajectory &traj, TAbs t)
    {
        auto segment_iter = upper_bound(traj._durations.begin(), traj._durations.end(), t);
        auto segment_idx = distance(traj._durations.begin(), segment_iter) - 1;
        segment_idx = min((long)traj.size() - 1L, max(0L, segment_idx)); // clip the segment index

        auto t0 = traj._durations[segment_idx], t1 = traj._durations[segment_idx+1];

        return PathPosition(segment_idx, (t - t0) / (t1 - t0));
    }

    // this method is analog to PathPosition::from_s(path, s_list)
    vector<PathPosition> PathPosition::from_t(const Trajectory &path, const std::vector<TAbs> &t_list)
    {
        vector<PathPosition> result;
        result.reserve(t_list.size());

        TRel cur_t = 0;
        size_t cur_idx = 1; // index of first point that has s larger than cur_s
        for (auto t : t_list)
            if (t < cur_t)
                result.push_back(PathPosition::from_t(path, t));
            else
            {
                while (t > path._durations[cur_idx] && cur_idx < path.size())
                    cur_idx++;

                auto t0 = path._durations[cur_idx-1], t1 = path._durations[cur_idx];
                result.push_back(PathPosition(cur_idx - 1, (t - t0) / (t1 - t0)));

                cur_t = t;
            }

        return result;
    }

    Vector2 Trajectory::solve_velocity(size_t segment_idx) const
    {
        assert (segment_idx >= 0 && segment_idx <= _line.size() - 2);

        TRel dt = _durations[segment_idx+1] - _durations[segment_idx];
        Vector2 p1(_line[segment_idx+1].get<0>(), _line[segment_idx+1].get<1>());
        Vector2 p2(_line[segment_idx].get<0>(), _line[segment_idx].get<1>());
        assert(dt != 0);
        return (p1 - p2) / dt;
    }

    Vector2 Trajectory::velocity_at(const PathPosition &pos, bool interpolate) const
    {
        if (!interpolate ||
            (pos.segment == 0 && pos.fraction < 0.5) || // assume constant speed at both ends
            (pos.segment == _line.size() - 2 && pos.fraction >= 0.5))
            return solve_velocity(pos.segment);
        else
        {
            // linear interpolate based on fraction (position)
            Vector2 vel1, vel2;
            TRel w1, w2;
            if (pos.fraction < 0.5)
            {
                vel1 = solve_velocity(pos.segment - 1);
                w1 = 0.5 - pos.fraction;
                vel2 = solve_velocity(pos.segment);
                w2 = 0.5 + pos.fraction;
            }
            else
            {
                vel1 = solve_velocity(pos.segment);
                w1 = 1.5 - pos.fraction;
                vel2 = solve_velocity(pos.segment + 1);
                w2 = pos.fraction - 0.5;
            }
            return vel1.array() * w1 + vel2.array() * w2;
        }
    }

    Vector2 Trajectory::solve_acceleration(size_t point_idx) const
    {
        assert (point_idx >= 1 && point_idx <= _line.size() - 2);

        TRel dt = (_durations[point_idx+1] - _durations[point_idx-1]) / 2;
        Vector2 vel1 = solve_velocity(point_idx - 1);
        Vector2 vel2 = solve_velocity(point_idx);
        assert(dt != 0);
        return (vel2.array() - vel1.array()) / dt;
    }

    Vector2 Trajectory::acceleration_at(const PathPosition &pos, bool interpolate) const
    {
        if (_line.size() <= 2) return Vector2(0, 0); // no acceleration for one segment

        if (pos.segment == 0)
            return solve_acceleration(1);
        else if (pos.segment == _line.size() - 2)
            return solve_acceleration(pos.segment);
        else if (interpolate)
        {
            return solve_acceleration(pos.segment).array() * (1 - pos.fraction) +
                   solve_acceleration(pos.segment + 1).array() * pos.fraction;
        }
        else
        {
            if (pos.fraction < 0.5)
                return solve_acceleration(pos.segment);
            else
                return solve_acceleration(pos.segment+1);
        }
    }

    Trajectory Trajectory::resample_at(const std::vector<TAbs> &t_list) const
    {
        vector<PathPosition> pos_list = PathPosition::from_t(*this, t_list);
        vector<Point> plist; plist.reserve(t_list.size());
        transform(pos_list.begin(), pos_list.end(),
            std::back_inserter(plist), std::bind(&Path::point_at, this, _1));
        return Trajectory(Path(move(plist)), t_list);
    }

    QuinticPolyTrajectory::QuinticPolyTrajectory(
        TRel T, const Vector3 &x0, const Vector3 &xT,
        const Vector3 &y0, const Vector3 &yT, bool relax_sx
    ) : _x_coeffs(), _y_coeffs(), _T(T)
    {
        assert (T > 0);
        TRel T3 = T * T * T;

        Vector3 c012x; c012x << x0(0), x0(1), x0(2) / 2;
        Vector3 c012y; c012y << y0(0), y0(1), y0(2) / 2;

        Matrix3 M1, M2;
        M1 << 1, T, T*T,
              0, 1, 2*T,
              0, 0, 2;
        M2 << T3, T3*T, T3*T*T,
              3*T*T, 4*T3, 5*T3*T,
              6*T, 12*T*T, 20*T3;
        Matrix3 M2inv = M2.inverse();

        auto c345y = M2inv * (yT - M1 * c012y);
        _y_coeffs << c345y(2), c345y(1), c345y(0), c012y(2), c012y(1), c012y(0);

        if (relax_sx)
        {
            auto c34x = M2.block(1,0,2,2).inverse() * (xT.tail(2) - M1.block(1,1,2,2) * c012x.tail(2));
            _x_coeffs << 0, c34x(1), c34x(0), c012x(2), c012x(1), c012x(0);
        }
        else
        {
            auto c345x = M2inv * (xT - M1 * c012x);
            _x_coeffs << c345x(2), c345x(1), c345x(0), c012x(2), c012x(1), c012x(0);
        }
    }

    Point QuinticPolyTrajectory::point_at(TRel t) const
    {
        TRel x = _x_coeffs(0), y = _y_coeffs(0);
        for (size_t i = 1; i < 6; i++)
        {
            x = x * t + _x_coeffs(i);
            y = y * t + _y_coeffs(i);
        }
        return Point(x, y);
    }

    Vector2 QuinticPolyTrajectory::velocity_at(TRel t) const
    {
        TRel x = 5 * _x_coeffs(0), y = 5 * _y_coeffs(0);
        for (size_t i = 1; i < 5; i++)
        {
            x = x * t + (5-i) * _x_coeffs(i);
            y = y * t + (5-i) * _y_coeffs(i);
        }
        return Vector2(x, y);
    }

    Vector2 QuinticPolyTrajectory::acceleration_at(TRel t) const
    {
        TRel x = 20 * _x_coeffs(0), y = 20 * _y_coeffs(0);
        for (size_t i = 1; i < 4; i++)
        {
            x = x * t + (5-i) * (4-i) * _x_coeffs(i);
            y = y * t + (5-i) * (4-i) * _y_coeffs(i);
        }
        return Vector2(x, y);
    }

    TRel QuinticPolyTrajectory::tangent_at(TRel t) const
    {
        auto vel = velocity_at(t);
        return atan2(vel(1), vel(0));
    }

    Trajectory QuinticPolyTrajectory::periodize(TRel interval) const
    {
        auto t = VectorX::LinSpaced((size_t)ceil(_T / interval) + 1, 0, _T).array();

        ArrayX x(t.rows()), y(t.rows());
        x.setConstant(_x_coeffs(0));
        y.setConstant(_y_coeffs(0));
        for (size_t i = 1; i < 6; i++)
        {
            x = x * t + _x_coeffs(i);
            y = y * t + _y_coeffs(i);
        }

        Trajectory result;
        result._line.reserve(t.rows());
        result._durations.reserve(t.rows());
        for (size_t i = 0; i < t.rows(); i++)
        {
            result._line.emplace_back(x(i), y(i));
            result._durations.emplace_back(t(i));
        }
        result.update_distance();
        return result;
    }

    Point CTRATrajectory::point_at(TRel t) const
    {
        TRel x = _init_state(0), y = _init_state(1), th = _init_state(2);
        TRel v = _init_state(3), a = _init_state(4), w = _init_state(5);

        TRel nth = th + w * t;
        TRel nv = v + a * t;
        TRel nx, ny;
        if (w == 0)
        {
            nx = x + (nv + v)/2 * cos(th) * t;
            ny = y + (nv + v)/2 * sin(th) * t;
        }
        else
        {
            nx = x + ( nv*w*sin(nth) + a*cos(nth) - v*w*sin(th) - a*cos(th)) / (w*w);
            ny = y + (-nv*w*cos(nth) + a*sin(nth) + v*w*cos(th) - a*sin(th)) / (w*w);
        }
        return Point(nx, ny);
    }

    Vector2 CTRATrajectory::velocity_at(TRel t) const
    {
        TRel nth = tangent_at(t);
        TRel nv = _init_state(3) + _init_state(4) * t;
        return Vector2(nv * cos(nth), nv * sin(nth));
    }

    Trajectory CTRATrajectory::periodize(TRel interval) const
    {
        auto t = VectorX::LinSpaced((size_t)ceil(_T / interval) + 1, 0, _T).array();

        TRel x = _init_state(0), y = _init_state(1), th = _init_state(2);
        TRel v = _init_state(3), a = _init_state(4), w = _init_state(5);

        auto nth = th + w * t;
        auto nv = v + a * t;
        VectorX nx, ny;
        if (w == 0)
        {
            nx = x + (nv + v)/2 * cos(th) * t;
            ny = y + (nv + v)/2 * sin(th) * t;
        }
        else
        {
            nx = x + ( nv*w*sin(nth) + a*cos(nth) - v*w*sin(th) - a*cos(th)) / (w*w);
            ny = y + (-nv*w*cos(nth) + a*sin(nth) + v*w*cos(th) - a*sin(th)) / (w*w);
        }

        Trajectory result;
        result._line.reserve(t.rows());
        result._durations.reserve(t.rows());
        for (size_t i = 0; i < t.rows(); i++)
        {
            result._line.emplace_back(nx(i), ny(i));
            result._durations.emplace_back(t(i));
        }
        result.update_distance();
        return result;
    }

    TRel tdistance(const Trajectory &lhs, const Trajectory &rhs)
    {
        // TODO: this implementation is not correct, we need 3D segment distance

        TRel min_dist = numeric_limits<TRel>::max();
        vector<PathPosition> rhs_pos = PathPosition::from_t(lhs, rhs.timestamps());
        for (int i = 0; i < rhs_pos.size(); i++)
        {
            TRel dist = distance(lhs.point_at(rhs_pos[i]), rhs.vertices()[i]);
            min_dist = min(min_dist, dist);
        }

        vector<PathPosition> lhs_pos = PathPosition::from_t(rhs, lhs.timestamps());
        for (int i = 0; i < lhs_pos.size(); i++)
        {
            TRel dist = distance(rhs.point_at(lhs_pos[i]), lhs.vertices()[i]);
            min_dist = min(min_dist, dist);
        }
        return min_dist;
    }
}

namespace std
{
    string to_string(const traji::Trajectory &value)
    {
        if (value.size() == 0)
            return string("[]");

        stringstream ss;
        ss << '[' << to_string(value.vertices()[0]) << " @ " << value.timestamps()[0];
        for (size_t i = 1; i < value.vertices().size(); i++)
            ss << ", " << to_string(value.vertices()[i]) << " @ " << value.timestamps()[i];
        ss << ']';
        return ss.str();
    }

    string to_string(const traji::QuinticPolyTrajectory &value)
    {
        stringstream ss;
        ss << "(T=" << value.T() << ", x_coeffs [" << value.x_coeffs()(0);
        for (size_t i = 1; i < 6; i++)
            ss << ", " << value.x_coeffs()(i);
        ss << "], y_coeffs [" << value.y_coeffs()(0);
        for (size_t i = 1; i < 6; i++)
            ss << ", " << value.y_coeffs()(i);
        ss << "])";
        return ss.str();
    }
} // namespace std
