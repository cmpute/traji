ctypedef float TFloat

cdef extern from "traji.hpp" namespace "traji":
    cdef cppclass Point:
        Point ()

    cdef cppclass Path:
        Path ()
        Point point_from(TFloat s)
