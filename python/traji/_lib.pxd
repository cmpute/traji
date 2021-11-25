from traji.decl cimport Path as cPath, Point as cPoint

cdef class Point:
    cdef cPoint _data

cdef class Path:
    cdef cPath _data
