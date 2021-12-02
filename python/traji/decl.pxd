from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.utility cimport pair

ctypedef float TFloat

cdef extern from "traji.hpp" namespace "traji":
    """
namespace traji {
    Vector6 make_vec6(TFloat v0, TFloat v1, TFloat v2,
                      TFloat v3, TFloat v4, TFloat v5) {
        Vector6 vec; vec << v0, v1, v2, v3, v4, v5; return vec;
    }
}
    """

    cdef cppclass Vector3:
        Vector3() except +
        Vector3(TFloat*) except + 
        Vector3(TFloat, TFloat, TFloat) except + 
        float *data()
        float& at "operator()"(size_t index)
        @staticmethod
        Vector3 Zero()

    cdef cppclass Vector6:
        Vector6() except +
        Vector6(TFloat*) except + 
        Vector6(TFloat, TFloat, TFloat) except + 
        float *data()
        float& at "operator()"(size_t index)
        @staticmethod
        Vector6 Zero()

    Vector6 make_vec6(TFloat,TFloat,TFloat,TFloat,TFloat,TFloat)

cdef extern from "traji.hpp" namespace "traji":
    cdef cppclass SegmentType:
        pass
    cdef SegmentType SegmentType_Line "traji::SegmentType::Line"
    cdef SegmentType SegmentType_Arc "traji::SegmentType::Arc"

    cdef cppclass Point:
        Point ()
        Point (TFloat x, TFloat y)
        TFloat x "get<0>" ()
        TFloat y "get<1>" ()

        bint operator==(Point)
        bint operator!=(Point)

    cdef cppclass PathPosition:
        PathPosition()
        size_t component
        size_t segment
        TFloat fraction

        bint operator==(PathPosition)
        bint operator!=(PathPosition)

        TFloat to_s(const Path &path)
        TFloat to_s(const HeteroPath &path)

        @staticmethod
        PathPosition from_s(const Path &path, TFloat s)
        @staticmethod
        PathPosition from_s_hetero "from_s" (const HeteroPath &path, TFloat s)
        @staticmethod
        vector[PathPosition] from_s_batch "from_s" (const Path &path, const vector[TFloat] &s)

        TFloat to_t(const Trajectory &traj)
        @staticmethod
        PathPosition from_t(const Trajectory &traj, TFloat t)

    cdef cppclass Path:
        Path ()
        Path (const Path&)
        Path (vector[Point].iterator begin, vector[Point].iterator end)

        size_t size()
        TFloat length()
        vector[TFloat] segment_lengths()

        vector[Point]& data()
        Point& operator[](size_t)

        bint operator==(Path)
        bint operator!=(Path)

        Point point_from(TFloat s)
        Point point_at(const PathPosition &pos)

        TFloat tangent_from(TFloat s)
        TFloat tangent_at(const PathPosition &pos)

        TFloat interpolate_from(const vector[TFloat] &values, TFloat s)
        TFloat interpolate_at(const vector[TFloat] &values, const PathPosition &pos)

        pair[TFloat, PathPosition] project(const Point &point)

        Path respacing(TFloat resolution, TFloat smooth_radius)
        Path densify(TFloat resolution)
        HeteroPath smooth(TFloat smooth_radius, SegmentType segtype)
        Path resample(const vector[TFloat] &s_list)

    cdef cppclass Trajectory(Path):
        Trajectory()
        Trajectory(const Trajectory& traj)
        Trajectory(const Path& path, vector[TFloat] &timestamps)
        Trajectory(vector[Point].iterator begin, vector[Point].iterator end,
                   vector[TFloat].iterator t_begin, vector[TFloat].iterator t_end)

        const vector[TFloat]& timestamps()

        Point point_at(TFloat t)

    cdef cppclass QuinticPolyTrajectory:
        QuinticPolyTrajectory(TFloat T, const Vector6 &x_coeffs, const Vector6 &y_coeffs)
        QuinticPolyTrajectory(TFloat T, const Vector3 &x0, const Vector3 &xT,
            const Vector3 &y0, const Vector3 &yT)

        const Vector6& x_coeffs()
        const Vector6& y_coeffs()
        TFloat T()

        Point point_at(TFloat t)
        TFloat tangent_at(TFloat t)

        Trajectory rasterize(TFloat resolution)
        Trajectory periodize(TFloat interval)

    cdef cppclass HeteroPath:
        HeteroPath()
        HeteroPath(const HeteroPath&)

cdef extern from "traji.hpp" namespace "traji::frenet":
    Point from_cartesian(const Path &ref, const Point &point)
    Point to_cartesian(const Path &ref, const Point &point)

    Path from_cartesian(const Path &ref, const Path &path)
    Path to_cartesian(const Path &ref, const Path &path)

    Trajectory from_cartesian(const Path &ref, const Trajectory &traj)
    Trajectory to_cartesian(const Path &ref, const Trajectory &traj)

cdef extern from "traji.hpp" namespace "std":
    string to_string(const Point &value)
    string to_string(const PathPosition &value)
    string to_string(const Path &value)
    string to_string(const Trajectory &value)
    string to_string(const QuinticPolyTrajectory &value)
