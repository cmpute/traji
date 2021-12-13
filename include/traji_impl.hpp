#pragma once

#include "traji.hpp"
using namespace std;

namespace traji
{
    /// Some helper functions for line segments
    namespace line
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
        inline pair<TFloat, TFloat> sdistance(const Point &p0, const Point &p1, const TFloat l, const Point &p)
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
        inline ArcParams solve_smooth(
            const Point &p0, const Point &pivot, const Point &p1,
            TFloat l0, TFloat l1, TFloat smooth_radius)
        {
            auto tan1 = line::tangent(p0, pivot);
            auto tan2 = line::tangent(pivot, p1);

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

        inline TFloat tangent(const ArcParams &params, TFloat fraction)
        {
            auto angular_pos = params.start_angle + params.angle * fraction;
            return angular_pos + (params.angle > 0 ? pi2 : -pi2);
        }
    }
}
