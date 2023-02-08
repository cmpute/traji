#include <stdexcept>
#include "traji.hpp"
#include "traji_impl.hpp"

namespace traji
{
    void HeteroPath::update_distance(TRel s0)
    {
        if (_points.size() == 1) // Fix invalid line string (only 1 point)
        {
            _points.clear();
            return;
        }

        TRel s = s0;
        _distance.resize(_points.size());
        _distance[0] = s;
        for (int i = 1; i < _points.size(); i++)
        {
            s += distance(_points[i], _points[i-1]);
            _distance[i] = s;
        }
    }

    arc::ArcParams to_arc (const HeteroSegment &segment,
        const Point &start, const Point &end)
    {
        Point center = Point(segment.params[0], segment.params[1]);
        TRel radius = distance(center, start);

        TRel tan1 = line::tangent(center, start);
        TRel tan2 = line::tangent(center, end);
        auto turn_angle = arc::warp_angle(tan2 - tan1);

        return arc::ArcParams { center, radius, tan2 - tan1, tan1 };
    }

    TRel HeteroSegment::length (const Point &start, const Point &end) const
    {
        switch (type)
        {
            case SegmentType::Arc:
            {
                auto params = to_arc(*this, start, end);
                return abs(params.angle) * params.radius;
            }
            case SegmentType::Line:
                return distance(start, end);
            default:
                throw invalid_argument("segment type");
        }
    }

    Point HeteroSegment::point_at (const Point &start, const Point &end, TRel fraction) const
    {
        switch (type)
        {
            case SegmentType::Arc:
            {
                auto params = to_arc(*this, start, end);
                return arc::interpolate(params, fraction);
            }
            case SegmentType::Line:
                return line::interpolate(start, end, fraction);
            default:
                throw invalid_argument("segment type");
        }
    }

    TRel HeteroSegment::tangent_at (const Point &start, const Point &end, TRel fraction) const
    {
        switch (type)
        {
            case SegmentType::Arc:
            {
                auto params = to_arc(*this, start, end);
                return arc::tangent(params, fraction);
            }
            case SegmentType::Line:
                return line::tangent(start, end);
            default:
                throw;
        }
    }

    TAbs PathPosition::to_s(const HeteroPath& path) const
    {
        return path._distance[segment] + (path._distance[segment+1] - path._distance[segment]) * fraction;
    }

    PathPosition PathPosition::from_s(const HeteroPath& path, TAbs s)
    {
        auto segment_iter = upper_bound(path._distance.begin(), path._distance.end(), s);
        auto segment_idx = distance(path._distance.begin(), segment_iter) - 1;
        segment_idx = min((long)path.size() - 1L, max(0L, segment_idx)); // clip the segment index

        auto s0 = path._distance[segment_idx], s1 = path._distance[segment_idx+1];

        return PathPosition(segment_idx, (s - s0) / (s1 - s0));
    }

    Point HeteroPath::point_at(const PathPosition &pos)
    {
        return _segments[pos.segment].point_at(
            _points[pos.segment], _points[pos.segment+1], pos.fraction);
    }

    TRel HeteroPath::tangent_at(const PathPosition &pos)
    {
        return _segments[pos.segment].tangent_at(
            _points[pos.segment], _points[pos.segment+1], pos.fraction);
    }

    Path HeteroPath::rasterize(TRel resolution) const // similar to Path::densify()
    {
        Path result;
        if (size() == 0)
            return result;

        result._line.push_back(_points.front());
        for (size_t i = 1; i < _points.size(); i++)
        {
            int mul = ceil((_distance[i] - _distance[i-1]) / resolution);
            if (mul <= 1)
                result._line.push_back(_points[i]);
            else
            {
                for (size_t j = 1; j <= mul; j++)
                    result._line.emplace_back(
                        _segments[i-1].point_at(_points[i-1], _points[i], (TRel) j / mul));
            }
        }

        result.update_distance();
        result.s0 = _distance[0];
        return result;
    }
}
