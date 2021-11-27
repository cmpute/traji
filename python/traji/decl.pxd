from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.utility cimport pair

ctypedef float TFloat

cdef extern from "traji.hpp" namespace "traji":
    cdef cppclass Point:
        Point ()
        Point (TFloat x, TFloat y)
        TFloat x "get<0>" ()
        TFloat y "get<1>" ()

    cdef cppclass PathPosition:
        PathPosition()
        size_t component
        size_t segment
        TFloat fraction

        TFloat to_s(const Path &path)

        @staticmethod
        PathPosition from_s(const Path &path, TFloat s)

    cdef cppclass Path:
        Path ()
        Path (vector[Point].iterator begin, vector[Point].iterator end)

        size_t size()
        TFloat length()
        vector[TFloat] segment_lengths()

        vector[Point]& data()
        Point& operator[](size_t)

        Point point_from(TFloat s)
        Point point_at(const PathPosition &pos)

        TFloat tangent_from(TFloat s)
        TFloat tangent_at(const PathPosition &pos)

        TFloat interpolate_from(const vector[TFloat] &values, TFloat s)
        TFloat interpolate_from(const vector[TFloat] &values, const PathPosition &pos)

        pair[TFloat, PathPosition] project(const Point &point)

        Path respacing(TFloat resolution, TFloat smooth_radius)
        Path densify(TFloat resolution)
        Path smooth(TFloat resolution, TFloat smooth_radius)

cdef extern from "traji.hpp" namespace "std":
    string to_string(const Point &value)
    string to_string(const PathPosition &value)
    string to_string(const Path &value)
