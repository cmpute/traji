
cdef class Point:
    property x:
        def __get__(self):
            return self._data.x()
    property y:
        def __get__(self):
            return self._data.y()

cdef class Path:
    def __cinit__(self):
        pass

