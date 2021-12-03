#include <cmath>
#include <Eigen/Dense>
#include "traji.hpp"

namespace traji
{
    TFloat PathPosition::to_t(const Trajectory &traj)
    {
        traj._timestamps[segment] + (traj._timestamps[segment+1] - traj._timestamps[segment]) * fraction;
    }
    PathPosition PathPosition::from_t(const Trajectory &traj, TFloat t)
    {
        auto segment_iter = lower_bound(traj._timestamps.begin(), traj._timestamps.end(), t);
        auto segment_idx = distance(traj._timestamps.begin(), segment_iter);
        auto t0 = traj._timestamps[segment_idx], t1 = traj._timestamps[segment_idx+1];

        PathPosition result;
        result.segment = segment_idx;
        result.fraction = (t - t0) / (t1 - t0);
        return result;
    }

    QuinticPolyTrajectory::QuinticPolyTrajectory(
        TFloat T, const Vector3 &x0, const Vector3 &xT,
        const Vector3 &y0, const Vector3 &yT
    ) : _x_coeffs(), _y_coeffs(), _T(T)
    {
        assert (T > 0);
        TFloat T3 = T * T * T;

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

        auto c345x = M2inv * (xT - M1 * c012x);
        auto c345y = M2inv * (yT - M1 * c012y);

        _x_coeffs << c345x(2), c345x(1), c345x(0), c012x(2), c012x(1), c012x(0);
        _y_coeffs << c345y(2), c345y(1), c345y(0), c012y(2), c012y(1), c012y(0);
    }

    Point QuinticPolyTrajectory::point_at(TFloat t) const
    {
        TFloat x = _x_coeffs(0), y = _y_coeffs(0);
        for (size_t i = 1; i < 6; i++)
        {
            x = x * t + _x_coeffs(i);
            y = y * t + _y_coeffs(i);
        }
        return Point(x, y);
    }

    Vector2 QuinticPolyTrajectory::velocity_at(TFloat t) const
    {
        TFloat x = 5 * _x_coeffs(0), y = 5 * _y_coeffs(0);
        for (size_t i = 1; i < 5; i++)
        {
            x = x * t + (5-i) * _x_coeffs(i);
            y = y * t + (5-i) * _y_coeffs(i);
        }
        return Vector2(x, y);
    }

    TFloat QuinticPolyTrajectory::tangent_at(TFloat t) const
    {
        auto vel = velocity_at(t);
        return atan2(vel(1), vel(0));
    }

    Trajectory QuinticPolyTrajectory::periodize(TFloat interval) const
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
        result._timestamps.reserve(t.rows());
        for (size_t i = 0; i < t.rows(); i++)
        {
            result._line.emplace_back(x(i), y(i));
            result._timestamps.emplace_back(t(i));
        }
        result.update_distance();
        return result;
    }
}

namespace std
{
    std::string to_string(const traji::Trajectory &value)
    {
        if (value.size() == 0)
            return string("[]");

        stringstream ss;
        ss << '[' << to_string(value[0]) << " @ " << value.timestamps()[0];
        for (size_t i = 1; i < value.size(); i++)
            ss << ", " << to_string(value[i]) << " @ " << value.timestamps()[i];
        ss << ']';
        return ss.str();
    }

    std::string to_string(const traji::QuinticPolyTrajectory &value)
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
