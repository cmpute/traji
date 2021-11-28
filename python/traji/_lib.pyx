from libcpp.vector cimport vector
from libcpp.utility cimport pair
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc

from traji.decl cimport to_string, TFloat

cdef class Point:
    def __init__(self, x, y=None, bint _noinit=False): # accept Point(x, y) or Point([x, y])
        if not _noinit:
            if isinstance(x, (int, float)) and isinstance(y, (int, float)):
                self._data = cPoint(x, y)
            elif isinstance(x, (tuple, list)):
                self._data = cPoint(x[0], x[1])
            else:
                raise ValueError("Unrecognized input!")

    @staticmethod
    cdef Point wrap(const cPoint &value):
        cdef Point p = Point(0, _noinit=True)
        p._data = value
        return p

    property x:
        def __get__(self): return self._data.x()
    property y:
        def __get__(self): return self._data.y()

    def __str__(self):
        return to_string(self._data).decode()
    def __repr__(self):
        return "<Point %s>" % to_string(self._data).decode()
    def __eq__(self, other):
        return isinstance(other, Point) and self._data == (<Point>other)._data
    def __ne__(self, other):
        return not (self == other)

cdef class PathPosition:
    def __init__(self, size_t component=0, size_t segment=0, TFloat fraction=0, bint _noinit=False):
        if not _noinit:
            self._data.component = component
            self._data.segment = segment
            self._data.fraction = fraction

    @staticmethod
    cdef PathPosition wrap(const cPathPosition &value):
        cdef PathPosition p = PathPosition(_noinit=True)
        p._data = value
        return p

    property component:
        def __get__(self): return self._data.component
        def __set__(self, size_t val): self._data.component = val
    property segment:
        def __get__(self): return self._data.segment
        def __set__(self, size_t val): self._data.segment = val
    property fraction:
        def __get__(self): return self._data.fraction
        def __set__(self, size_t val): self._data.fraction = val

    def __str__(self):
        return to_string(self._data).decode()
    def __repr__(self):
        return "<PathPosition %d, %d, %.3f>" % (self._data.component, self._data.segment, self._data.fraction)
    def __eq__(self, other):
        return isinstance(other, PathPosition) and self._data == (<PathPosition>other)._data
    def __ne__(self, other):
        return not (self == other)

cdef class _PathIterator:
    cdef vector[cPoint].iterator it, end

    def __next__(self):
        if self.it == self.end:
            raise StopIteration
        else:
            p = Point.wrap(deref(self.it))
            preinc(self.it)
            return p

cdef class Path:
    def __init__(self, points, bint _noinit=False):
        cdef vector[cPoint] cpoints
        if not _noinit and points:
            cpoints = vector[cPoint](len(points))
            for i, p in enumerate(points):
                if isinstance(p, Point):
                    cpoints[i] = (<Point>p)._data
                else:
                    cpoints[i] = cPoint(p[0], p[1])
            self._data = cPath(cpoints.begin(), cpoints.end())

    @staticmethod
    cdef Path wrap(const cPath &value):
        cdef Path p = Path(None, _noinit=True)
        p._data = value
        return p

    def __len__(self):
        return self._data.size()
    def __str__(self):
        return to_string(self._data).decode()
    def __repr__(self):
        return "<Path with %d points>" % self._data.size()
    def __eq__(self, other):
        return isinstance(other, Path) and self._data == (<Path>other)._data
    def __ne__(self, other):
        return not (self == other)
    def __iter__(self):
        cdef _PathIterator it = _PathIterator()
        it.it = self._data.data().begin()
        it.end = self._data.data().end()
        return it
    def __getitem__(self, size_t index):
        cdef Point p = Point(0, _noinit=True)
        p._data = self._data.data()[index]
        return p

    property length:
        def __get__(self): return self._data.length()
    property segment_lengths:
        def __get__(self): return self._data.segment_lengths()

    def point_from(self, TFloat s):
        return Point.wrap(self._data.point_from(s))
    def point_at(self, PathPosition pos):
        return Point.wrap(self._data.point_at(pos._data))

    def tangent_from(self, TFloat s):
        return self._data.tangent_from(s)
    def tangent_at(self, PathPosition pos):
        return self._data.tangent_at(pos._data)

    def interpolate_from(self, values, TFloat s):
        pass
    def interpolate_at(self, values, PathPosition pos):
        pass

    def project(self, Point point):
        cdef pair[TFloat, cPathPosition] result = self._data.project(point._data)
        return result.first, PathPosition.wrap(result.second)

    def respacing(self, TFloat resolution, TFloat smooth_radius=0):
        return Path.wrap(self._data.respacing(resolution, smooth_radius))
    def densify(self, TFloat resolution):
        return Path.wrap(self._data.densify(resolution))
    def smooth(self, TFloat resolution, TFloat smooth_radius=0):
        return Path.wrap(self._data.smooth(resolution, smooth_radius))
