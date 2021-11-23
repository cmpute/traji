from traji cimport Path as cPath

cdef class Path:
    cdef cPath data
    def __cinit__(self):
        pass

