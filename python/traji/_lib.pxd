from traji.decl cimport Path as cPath, Point as cPoint, PathPosition as cPathPosition

cdef class Point:
    cdef cPoint _data

    @staticmethod
    cdef Point wrap(const cPoint &value)

cdef class PathPosition:
    cdef cPathPosition _data
    
    @staticmethod
    cdef PathPosition wrap(const cPathPosition &value)

cdef class Path:
    cdef cPath _data

    @staticmethod
    cdef Path wrap(const cPath &value)
