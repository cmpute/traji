#include "traji.hpp"

namespace traji
{
    TFloat HeteroSegment::length (const Point &start, const Point &end) const
    {
        switch (type)
        {
            case SegmentType::Arc:
                break;
            case SegmentType::Line:
                break;
            default:
                throw;
        }
    }

    Point HeteroSegment::point_at (const Point &start, const Point &end, TFloat fraction) const
    {
        switch (type)
        {
            case SegmentType::Arc:
                break;
            case SegmentType::Line:
                break;
            default:
                throw;
        }
    }

    Point HeteroSegment::tangent_at (const Point &start, const Point &end, TFloat fraction) const
    {
        switch (type)
        {
            case SegmentType::Arc:
                break;
            case SegmentType::Line:
                break;
            default:
                throw;
        }
    }

    TFloat PathPosition::to_s(const HeteroPath& path)
    {
        return path._distance[segment] + (path._distance[segment+1] - path._distance[segment]) * fraction;
    }

    PathPosition PathPosition::from_s(const HeteroPath& path, TFloat s)
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

    TFloat HeteroPath::tangent_at(const PathPosition &pos)
    {
        return _segments[pos.segment].tangent_at(
            _points[pos.segment], _points[pos.segment+1], pos.fraction);
    }

    Path HeteroPath::rasterize(TFloat resolution) const
    {
        
    }
}
