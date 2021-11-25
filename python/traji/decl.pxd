ctypedef float TFloat

cdef extern from "traji.hpp" namespace "traji":
    cdef cppclass Point:
        Point ()
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
        Point point_from(TFloat s)
        Point point_at(const PathPosition &pos)
