from libc.stdlib cimport malloc, free
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

    @classmethod
    def from_s(cls, path, s):
        cdef vector[TFloat] s_list
        cdef vector[cPathPosition] p_list
        if isinstance(path, Path):
            if isinstance(s, (int, float)):
                return PathPosition.wrap(cPathPosition.from_s(deref((<Path>path)._ptr), <TFloat>s))
            elif isinstance(s, (list, tuple)):
                s_list = s
                p_list = cPathPosition.from_s_batch(deref((<Path>path)._ptr), s_list)
                return [PathPosition.wrap(p) for p in p_list]
            else:
                raise ValueError("Invalid s input!")
        else:
            raise ValueError("Unsupported path type!")

    def to_s(self, path):
        if isinstance(path, Path):
            return self._data.to_s(deref((<Path>path)._ptr))
        else:
            raise ValueError("Unsupported path type!")

    def to_t(self, traj):
        if isinstance(traj, Trajectory):
            return self._data.to_s(deref((<Trajectory>traj)._ptr))
        else:
            raise ValueError("Unsupported path input!")

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
    def __cinit__(self, points, bint _noinit=False):
        cdef vector[cPoint] cpoints
        if type(self) is Path:
            if _noinit:
                self._ptr = NULL
            elif points:
                cpoints = vector[cPoint](len(points))
                for i, p in enumerate(points):
                    if isinstance(p, Point):
                        cpoints[i] = (<Point>p)._data
                    else:
                        cpoints[i] = cPoint(p[0], p[1])
                self._ptr = new cPath(cpoints.begin(), cpoints.end())
            else:
                self._ptr = new cPath()

    def __dealloc__(self):
        if type(self) is Path:
            del self._ptr

    @staticmethod
    cdef Path wrap(const cPath &value):
        cdef Path p = Path(None, _noinit=True)
        p._ptr = new cPath(value)
        return p

    def __len__(self):
        return self._ptr.size()
    def __str__(self):
        return to_string(deref(self._ptr)).decode()
    def __repr__(self):
        return "<Path with %d points>" % self._ptr.size()
    def __eq__(self, other):
        return isinstance(other, Path) and self._ptr == (<Path>other)._ptr
    def __ne__(self, other):
        return not (self == other)
    def __iter__(self):
        cdef _PathIterator it = _PathIterator()
        it.it = self._ptr.data().begin()
        it.end = self._ptr.data().end()
        return it
    def __getitem__(self, size_t index):
        cdef Point p = Point(0, _noinit=True)
        p._data = self._ptr.data()[index]
        return p

    property length:
        def __get__(self): return self._ptr.length()
    property segment_lengths:
        def __get__(self): return self._ptr.segment_lengths()

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef Py_ssize_t *shape = <Py_ssize_t*>malloc(2*sizeof(Py_ssize_t))
        cdef Py_ssize_t *strides = <Py_ssize_t*>malloc(2*sizeof(Py_ssize_t))
        shape[0] = self._ptr.size()
        shape[1] = 2
        strides[0] = sizeof(TFloat) * 2
        strides[1] = sizeof(TFloat)

        buffer.buf = self._ptr.data().data()
        buffer.format = 'f' if sizeof(TFloat) == 4 else 'd'
        buffer.internal = NULL
        buffer.itemsize = sizeof(TFloat)
        buffer.len = sizeof(TFloat) * 2
        buffer.ndim = 2
        buffer.obj = self
        buffer.readonly = 1
        buffer.shape = shape
        buffer.strides = strides
        buffer.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buffer):
        free(buffer.shape)
        free(buffer.strides)

    def point_from(self, TFloat s):
        return Point.wrap(self._ptr.point_from(s))
    def point_at(self, PathPosition pos):
        return Point.wrap(self._ptr.point_at(pos._data))

    def tangent_from(self, TFloat s):
        return self._ptr.tangent_from(s)
    def tangent_at(self, PathPosition pos):
        return self._ptr.tangent_at(pos._data)

    def interpolate_from(self, values, TFloat s):
        pass
    def interpolate_at(self, values, PathPosition pos):
        pass

    def project(self, Point point):
        cdef pair[TFloat, cPathPosition] result = self._ptr.project(point._data)
        return result.first, PathPosition.wrap(result.second)

    def respacing(self, TFloat resolution, TFloat smooth_radius=0):
        return Path.wrap(self._ptr.respacing(resolution, smooth_radius))
    def densify(self, TFloat resolution):
        return Path.wrap(self._ptr.densify(resolution))

cdef class Trajectory:
    def __cinit__(self, points, timestamps=None, bint _noinit=False):
        cdef vector[cPoint] cpoints
        cdef vector[TFloat] tstamps

        if type(self) is Trajectory:
            if _noinit:
                self._ptr = NULL
            elif isinstance(points, Path):
                if not timestamps:
                    raise ValueError("Timestamps are required")
                tstamps = timestamps
                self._ptr = new cTrajectory(deref((<Path>points)._ptr), tstamps)
            elif points:
                cpoints = vector[cPoint](len(points))
                for i, p in enumerate(points):
                    if isinstance(p, Point):
                        cpoints[i] = (<Point>p)._data
                    else:
                        cpoints[i] = cPoint(p[0], p[1])
                if not isinstance(points[0], Point) and len(points[0]) == 3:
                    tstamps = vector[TFloat](len(points))
                    for i, p in enumerate(points):
                        tstamps[i] = p[2]
                if timestamps is not None:
                    tstamps = timestamps
                if tstamps.size() == 0:
                    raise ValueError("Timestamps are required")
                self._ptr = new cTrajectory(cpoints.begin(), cpoints.end(),
                    tstamps.begin(), tstamps.end())
            else:
                self._ptr = new cTrajectory()

    def __dealloc__(self):
        # only free the pointer at base class level
        if type(self) is Path:
            del self._ptr

    @staticmethod
    cdef Trajectory wrap(const cTrajectory &value):
        cdef Trajectory p = Trajectory(None, _noinit=True)
        p._ptr = new cTrajectory(value)
        return p
