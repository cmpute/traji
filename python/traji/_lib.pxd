from traji.decl cimport Path as cPath, Point as cPoint, \
    PathPosition as cPathPosition, Trajectory as cTrajectory

cdef class Point:
    cdef cPoint _data
    @staticmethod
    cdef Point wrap(const cPoint &value)

cdef class PathPosition:
    cdef cPathPosition _data
    @staticmethod
    cdef PathPosition wrap(const cPathPosition &value)

cdef class Path:
    cdef cPath* _ptr
    @staticmethod
    cdef Path wrap(const cPath &value)

cdef class Trajectory:
    cdef cTrajectory* _ptr
    @staticmethod
    cdef Trajectory wrap(const cTrajectory &value)
