# cython: profile=True
# TODO: the flag above is only for debug

from libc.stdlib cimport malloc, free
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc

from traji.decl cimport (to_string, TFloat, Vector2, Vector3, Vector6,
    make_vec6, SegmentType, SegmentType_Arc, SegmentType_Line,
    from_cartesian, to_cartesian,
    distance as cdistance, arg_distance as carg_distance,
    intersection as cintersection, arg_intersection as carg_intersection)

cdef class Point:
    def __init__(self, x, y=None, bint _noinit=False): # accept Point(x, y) or Point([x, y])
        if not _noinit:
            if isinstance(x, (int, float)) and isinstance(y, (int, float)):
                self._data = cPoint(x, y)
            elif isinstance(x, (tuple, list)):
                self._data = cPoint(x[0], x[1])
            elif hasattr(x, 'x') and hasattr(x, 'y'):
                self._data = cPoint(x.x, x.y)
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
    def __len__(self):
        return 2
    def __eq__(self, other):
        return isinstance(other, Point) and self._data == (<Point>other)._data
    def __ne__(self, other):
        return not (self == other)
    def __iter__(self):
        return iter((self._data.x(), self._data.y()))
    def __array__(self):
        return self.numpy()
    def numpy(self):
        import numpy as np
        return np.array([self._data.x(), self._data.y()])
    def shapely(self):
        from shapely.geometry import Point
        return Point(self._data.x(), self._data.y())

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

    @classmethod
    def from_t(cls, traj, t):
        if isinstance(traj, Trajectory):
            if isinstance(t, (int, float)):
                return PathPosition.wrap(cPathPosition.from_t(deref((<Trajectory>traj).ptr()), <TFloat>t))
            else:
                raise ValueError("Invalid t input!")
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
                if hasattr(points, 'coords'): # support for shapely geometries
                    points = points.coords

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
        return isinstance(other, Path) and deref(self._ptr) == deref((<Path>other)._ptr)
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

    def numpy(self):
        import numpy as np
        return np.asarray(self)

    def shapely(self):
        from shapely.geometry import LineString
        plist = [(p.x(), p.y()) for p in self._ptr.data()]
        return LineString(plist)

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
    def smooth(self, TFloat smooth_radius, str segtype = "arc"):
        cdef SegmentType ctype

        segtype = segtype.lower()
        if segtype == 'line':
            ctype = SegmentType_Line
        elif segtype == 'arc':
            ctype = SegmentType_Arc
        else:
            raise ValueError("The valid segment types are {line, arc}")

        return HeteroPath.wrap(self._ptr.smooth(smooth_radius, ctype))

cdef class Trajectory(Path):
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

    cdef cTrajectory* ptr(self):
        return <cTrajectory*>(self._ptr)

    @staticmethod
    cdef Trajectory wrap(const cTrajectory &value):
        cdef Trajectory p = Trajectory(None, _noinit=True)
        p._ptr = new cTrajectory(value)
        return p

    def __str__(self):
        return to_string(deref(self.ptr())).decode()
    def __repr__(self):
        return "<Trajectory with %d points>" % self.ptr().size()

    property timestamps:
        def __get__(self): return self.ptr().timestamps()

    def numpy(self):
        import numpy as np
        points = super().numpy()
        ret = np.empty((len(self), 3), dtype=points.dtype)
        ret[:, :2] = points
        ret[:, 2] = self.timestamps
        return ret

    def point_at(self, pos):
        if isinstance(pos, PathPosition):
            super().point_at(pos)
        elif isinstance(pos, float):
            return Point.wrap(self.ptr().point_at(<TFloat>pos))
        else:
            raise ValueError("Unrecognized position input!")

    def tangent_at(self, pos):
        if isinstance(pos, PathPosition):
            super().tangent_at(pos)
        elif isinstance(pos, float):
            return self.ptr().tangent_at(<TFloat>pos)
        else:
            raise ValueError("Unrecognized position input!")

    def velocity_from(self, float s):
        cdef Vector2 vel = self.ptr().velocity_from(s)
        return vel.at(0), vel.at(1)
    def velocity_at(self, pos):
        cdef Vector2 vel
        if isinstance(pos, PathPosition):
            vel = self.ptr().velocity_at((<PathPosition>pos)._data)
        elif isinstance(pos, float):
            vel = self.ptr().velocity_at(<TFloat>pos)
        else:
            raise ValueError("Unrecognized position input!")
        return vel.at(0), vel.at(1)

cdef class QuinticPolyTrajectory:
    def __cinit__(self, TFloat T, x0=None, xT=None, y0=None, yT=None, x_coeffs=None, y_coeffs=None):
        cdef Vector6 cx_coeffs, cy_coeffs
        cdef Vector3 cx0, cxT, cy0, cyT

        if not (x_coeffs is None) and not (y_coeffs is None):
            assert len(x_coeffs) >= 6 and len(y_coeffs) >= 6
            cx_coeffs = make_vec6(x_coeffs[0], x_coeffs[1], x_coeffs[2], x_coeffs[3], x_coeffs[4], x_coeffs[5])
            cy_coeffs = make_vec6(y_coeffs[0], y_coeffs[1], y_coeffs[2], y_coeffs[3], y_coeffs[4], y_coeffs[5])
            self._ptr = new cQuinticPolyTrajectory(T, cx_coeffs, cy_coeffs)
        elif not (x0 is None) and not (xT is None):
            assert len(x0) >= 3 and len(xT) >= 3 and len(y0) >= 3 and len(yT) >= 3
            cx0 = Vector3(x0[0], x0[1], x0[2])
            cxT = Vector3(xT[0], xT[1], xT[2])
            cy0 = Vector3(y0[0], y0[1], y0[2])
            cyT = Vector3(yT[0], yT[1], yT[2])
            self._ptr = new cQuinticPolyTrajectory(T, cx0, cxT, cy0, cyT)
        else:
            raise ValueError("No input is given for the QuinticPolyTrajectory")

    def __dealloc__(self):
        del self._ptr

    def __str__(self):
        return to_string(deref(self._ptr)).decode()
    def __repr__(self):
        return "<QuinticPolyTrajectory with T=%.4f>" % (self._ptr.T())

    property T:
        def __get__(self): return self._ptr.T()
    property x_coeffs:
        def __get__(self):
            return [self._ptr.x_coeffs().at(i) for i in range(6)]
    property y_coeffs:
        def __get__(self):
            return [self._ptr.y_coeffs().at(i) for i in range(6)]

    def point_at(self, TFloat t):
        return Point.wrap(self._ptr.point_at(t))
    def tangent_at(self, TFloat t):
        return self._ptr.tangent_at(t)
    def velocity_at(self, TFloat t):
        cdef Vector2 vel = self._ptr.velocity_at(t)
        return vel.at(0), vel.at(1)

    def acceleration_at(self, TFloat t):
        cdef Vector2 accel = self._ptr.acceleration_at(t)
        return accel.at(0), accel.at(1)

    def periodize(self, TFloat interval):
        return Trajectory.wrap(self._ptr.periodize(interval))

cdef class HeteroPath:
    def __cinit__(self, bint _noinit = False):
        if not _noinit:
            self._ptr = new cHeteroPath()

    def __dealloc__(self):
        del self._ptr

    @staticmethod
    cdef HeteroPath wrap(const cHeteroPath &value):
        cdef HeteroPath p = HeteroPath(_noinit=True)
        p._ptr = new cHeteroPath(value)
        return p

def frenet_from_cartesian(Path ref, target):
    if isinstance(target, Point):
        return Point.wrap(from_cartesian(deref(ref._ptr), (<Point>target)._data))
    elif isinstance(target, Trajectory):
        return Trajectory.wrap(from_cartesian(deref(ref._ptr), deref((<Trajectory>target).ptr())))
    elif isinstance(target, Path):
        return Path.wrap(from_cartesian(deref(ref._ptr), deref((<Path>target)._ptr)))
    else:
        raise ValueError("Invalid target type!")

def frenet_to_cartesian(Path ref, target):
    if isinstance(target, Point):
        return Point.wrap(to_cartesian(deref(ref._ptr), (<Point>target)._data))
    elif isinstance(target, Trajectory):
        return Trajectory.wrap(to_cartesian(deref(ref._ptr), deref((<Trajectory>target).ptr())))
    elif isinstance(target, Path):
        return Path.wrap(to_cartesian(deref(ref._ptr), deref((<Path>target)._ptr)))
    else:
        raise ValueError("Invalid target type!")

def distance(lhs, rhs):
    if isinstance(lhs, Path):
        if isinstance(rhs, Point):
            return cdistance(deref((<Path>lhs)._ptr),  (<Point>rhs)._data)
        else:
            raise ValueError("Unsupported type")
    elif isinstance(lhs, Point):
        if isinstance(rhs, Point):
            return cdistance((<Point>lhs)._data, (<Point>rhs)._data)
        elif isinstance(rhs, Path):
            return cdistance(deref((<Path>rhs)._ptr), (<Point>lhs)._data)
        else:
            raise ValueError("Unsupported type")
    else:
        raise ValueError("Unsupported type")

def intersection(lhs, rhs):
    cdef vector[cPoint] points_result
    if isinstance(lhs, Path):
        if isinstance(rhs, Path):
            points_result = cintersection(deref((<Path>lhs)._ptr), deref((<Path>rhs)._ptr))
            return [Point.wrap(p) for p in points_result]
        else:
            raise ValueError("Unsupported type")

def arg_intersection(lhs, rhs):
    cdef vector[pair[cPathPosition, cPathPosition]] points_result
    if isinstance(lhs, Path):
        if isinstance(rhs, Path):
            points_result = carg_intersection(deref((<Path>lhs)._ptr), deref((<Path>rhs)._ptr))
            return [(PathPosition.wrap(p.first), PathPosition.wrap(p.second)) for p in points_result]
        else:
            raise ValueError("Unsupported type")
