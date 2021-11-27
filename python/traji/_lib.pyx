from libcpp.vector cimport vector
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc

from traji.decl cimport to_string

cdef class Point:
    def __init__(self, x, y=None, bint _noinit=False): # accept Point(x, y) or Point([x, y])
        if not _noinit:
            if isinstance(x, float) and isinstance(y, float):
                self._data = cPoint(x, y)
            elif isinstance(x, (tuple, list)):
                self._data = cPoint(x[0], x[1])

    property x:
        def __get__(self):
            return self._data.x()
    property y:
        def __get__(self):
            return self._data.y()

    def __str__(self):
        return to_string(self._data).decode()

    def __repr__(self):
        return "<Point %s>" % to_string(self._data).decode()

cdef class _PathIterator:
    cdef vector[cPoint].iterator it, end

    def __next__(self):
        if self.it == self.end:
            raise StopIteration
        else:
            p = Point(0, _noinit=True)
            p._data = deref(self.it)
            preinc(self.it)
            return p

cdef class Path:
    def __init__(self, points):
        cdef vector[cPoint] cpoints = vector[cPoint](len(points)) 
        if not points:
            self._data = cPath()
        else:
            for i, p in enumerate(points):
                if isinstance(p, Point):
                    cpoints[i] = (<Point>p)._data
                else:
                    cpoints[i] = cPoint(p[0], p[1])
            self._data = cPath(cpoints.begin(), cpoints.end())

    def __len__(self):
        return self._data.size()

    def __str__(self):
        return to_string(self._data).decode()

    def __repr__(self):
        return "<Path with %d points>" % self._data.size()

    def __iter__(self):
        cdef _PathIterator it = _PathIterator()
        it.it = self._data.data().begin()
        it.end = self._data.data().end()
        return it

    def __getitem__(self, size_t index):
        cdef Point p = Point(0, _noinit=True)
        p._data = self._data.data()[index]
        return p
