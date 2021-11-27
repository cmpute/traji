from traji.decl cimport Path as cPath, Point as cPoint, PathPosition as cPathPosition

cdef class Point:
    cdef cPoint _data

cdef class PathPosition:
    cdef cPathPosition _cPathPosition

cdef class Path:
    cdef cPath _data
