#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include "traji.hpp"

namespace py = pybind11;
using namespace traji;
using namespace std;

inline Point cast_point(py::object point)
{
    if (py::hasattr(point, "x") && py::hasattr(point, "y"))
        return Point(point.attr("x").cast<TFloat>(), point.attr("y").cast<TFloat>());
    else
    {
        vector<TFloat> data = point.cast<vector<TFloat>>();
        if (data.size() != 2)
            throw py::value_error("The length of point values should be 2!");
        return Point(data[0], data[1]);
    }
}

PYBIND11_MODULE(_bindings, m) {

    py::enum_<SegmentType>(m, "SegmentType")
        .value("Line", SegmentType::Line)
        .value("Arc", SegmentType::Arc)
        .value("QuadraticBezier", SegmentType::QuadraticBezier)
        .value("CubicBezier", SegmentType::CubicBezier)
        .value("Polynomial", SegmentType::Polynomial)
        ;

    py::class_<Point>(m, "Point")
        .def(py::init<TFloat, TFloat>())
        .def(py::init<Point>())
        .def(py::init(&cast_point))
        .def_property_readonly("x", [](const Point &p) { return p.get<0>(); })
        .def_property_readonly("y", [](const Point &p) { return p.get<1>(); })
        .def("__str__", py::overload_cast<const Point&>(&to_string))
        .def("__repr__", [](const Point &p) { 
            stringstream ss; ss << "<Point " << to_string(p) << ">"; return ss.str();
        })
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("__array__", [](const Point &p) {
            py::array_t<TFloat> result({2});
            auto r = result.mutable_unchecked<1>();
            r(0) = p.get<0>(); r(1) = p.get<1>();
            return result;
        })
        .def("numpy", [](py::object p) { return p.attr("__array__")(); })
        .def("shapely", [](const Point &p) {
            py::object shapely = py::module::import("shapely.geometry");
            py::object shapely_point = shapely.attr("Point");
            return shapely_point(p.get<0>(), p.get<1>());
        })
        ;

    py::class_<PathPosition>(m, "PathPosition")
        .def(py::init<>())
        .def(py::init<size_t, TFloat>(), py::arg("segment"), py::arg("fraction"))
        .def(py::init<size_t, size_t, TFloat>(), py::arg("component"), py::arg("segment"), py::arg("fraction"))
        .def("__str__", py::overload_cast<const PathPosition&>(to_string))
        .def("__repr__", [](const PathPosition &pos) { 
            stringstream ss; ss << "<PathPosition ";
            ss << pos.component << ", ";
            ss << pos.segment << ", ";
            ss << pos.fraction << ">";
            return ss.str();
        })
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self < py::self)
        .def(py::self > py::self)
        .def(py::self <= py::self)
        .def(py::self >= py::self)
        .def_readwrite("component", &PathPosition::component)
        .def_readwrite("segment", &PathPosition::segment)
        .def_readwrite("fraction", &PathPosition::fraction)
        .def("to_s", py::overload_cast<const Path&>(&PathPosition::to_s, py::const_))
        .def("to_s", py::overload_cast<const HeteroPath&>(&PathPosition::to_s, py::const_))
        .def_static("from_s", py::overload_cast<const Path&, TFloat>(PathPosition::from_s))
        .def_static("from_s", py::overload_cast<const HeteroPath&, TFloat>(PathPosition::from_s))
        .def_static("from_s", py::overload_cast<const Path&, const vector<TFloat>&>(PathPosition::from_s))
        .def("to_t", &PathPosition::to_t)
        .def_static("from_t", py::overload_cast<const Trajectory&, TFloat>(PathPosition::from_t))
        .def_static("from_t", py::overload_cast<const Trajectory&, const vector<TFloat>&>(PathPosition::from_t))
        .def("forward", &PathPosition::forward)
        .def("backward", &PathPosition::backward)
        ;

    py::class_<Path>(m, "Path", py::buffer_protocol())
        .def(py::init<>())
        .def(py::init<const Path&>(), py::arg("path"))
        .def(py::init<const vector<Point>&>(), py::arg("vertices"))
        .def(py::init<const vector<Point>&, TFloat>(), py::arg("vertices"), py::arg("s0"))
        .def(py::init([](py::object points) {
            if (py::hasattr(points, "coords"))
                points = points.attr("coords");

            auto pylist = points.cast<vector<py::object>>();
            vector<Point> clist; clist.reserve(pylist.size());
            for (auto p : pylist)
                clist.push_back(cast_point(p));
            return make_unique<Path>(move(clist));
        }), py::arg("vertices"))
        .def("__len__", &Path::size)
        .def("__str__", py::overload_cast<const Path&>(&to_string))
        .def("__repr__", [](const Path& p) {
            stringstream ss; ss << "<Path with " << p.size() << " segments>"; return ss.str();
        })
        .def("__iter__", [](const Path& p) {
            return py::make_iterator(p.vertices().begin(), p.vertices().end());
        }, py::keep_alive<0, 1>())
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def_property_readonly("vertices", py::overload_cast<>(&Path::vertices, py::const_))
        .def_property_readonly("length", &Path::length)
        .def_property_readonly("segment_lengths", &Path::segment_lengths)
        .def_buffer([](Path& p) {
            return py::buffer_info(
                static_cast<void*>(p.data().data()),
                sizeof(TFloat), py::format_descriptor<TFloat>::format(),
                2, {p.size() + 1, size_t(2)}, {2 * sizeof(TFloat), sizeof(TFloat)}
#if PYBIND11_VERSION_MAJOR >= 2 && PYBIND11_VERSION_MAJOR >= 4 && PYBIND11_VERSION_MAJOR >= 4
                , /* readonly = */ true
#endif
            );
        })
        .def("numpy", [](py::object self) { return py::array(self); })
        .def("shapely", [](py::object self) {
            py::object shapely = py::module::import("shapely.geometry");
            py::object shapely_linestring = shapely.attr("LineString");
            return shapely_linestring(py::array(self));
        })
        .def("point_from", &Path::point_from)
        .def("point_at", &Path::point_at)
        .def("tangent_from", &Path::tangent_from)
        .def("tangent_at", &Path::tangent_at)
        .def("curvature_from", &Path::curvature_from)
        .def("curvature_at", &Path::curvature_at)
        .def("interpolate_from", &Path::interpolate_from)
        .def("interpolate_at", &Path::interpolate_at)
        .def("project", &Path::project)
        .def("respacing", &Path::respacing, py::arg("resolution"), py::arg("smooth_radius") = 0)
        .def("densify", &Path::densify, py::arg("resolution"))
        .def("smooth", &Path::smooth, py::arg("smooth_radius"), py::arg("seg_type") = SegmentType::Arc)
        .def("resample_from", &Path::resample_from)
        .def("extract_from", &Path::extract_from)
        ;

    py::class_<Trajectory, Path>(m, "Trajectory")
        .def(py::init<>())
        .def(py::init<const Trajectory&>())
        .def(py::init<const Path&, const vector<TFloat>&>())
        .def(py::init<const Path&, TFloat, TFloat>())
        .def("__str__", py::overload_cast<const Trajectory&>(&to_string))
        .def("__repr__", [](const Trajectory& t) {
            stringstream ss; ss << "<Trajectory with " << t.size() << " points>"; return ss.str();
        })
        .def_property_readonly("timestamps", &Trajectory::timestamps)
        .def("numpy", [](const Trajectory& t) {
            py::array_t<TFloat> result({t.size() + 1, size_t(3)});
            auto r = result.mutable_unchecked<2>();
            for (size_t i = 0; i <= t.size(); i++)
            {
                r(i, 0) = t.vertices()[i].get<0>();
                r(i, 1) = t.vertices()[i].get<1>();
                r(i, 2) = t.timestamps()[i];
            }
            return result;
        })
        .def("point_at", py::overload_cast<TFloat>(&Trajectory::point_at, py::const_))
        .def("point_at", py::overload_cast<const PathPosition &>(&Trajectory::point_at, py::const_))
        .def("tangent_at", py::overload_cast<TFloat>(&Trajectory::tangent_at, py::const_))
        .def("tangent_at", py::overload_cast<const PathPosition &>(&Trajectory::tangent_at, py::const_))
        .def("velocity_from", &Trajectory::velocity_from, py::arg("s"), py::arg("interpolate") = false)
        .def("velocity_at", py::overload_cast<TFloat, bool>(&Trajectory::velocity_at, py::const_), py::arg("t"), py::arg("interpolate") = false)
        .def("velocity_at", py::overload_cast<const PathPosition&, bool>(&Trajectory::velocity_at, py::const_), py::arg("pos"), py::arg("interpolate") = false)
        .def("acceleration_from", &Trajectory::acceleration_from, py::arg("s"), py::arg("interpolate") = false)
        .def("acceleration_at", py::overload_cast<TFloat, bool>(&Trajectory::acceleration_at, py::const_), py::arg("t"), py::arg("interpolate") = false)
        .def("acceleration_at", py::overload_cast<const PathPosition&, bool>(&Trajectory::acceleration_at, py::const_), py::arg("pos"), py::arg("interpolate") = false)
        .def("resample_at", &Trajectory::resample_at)
        ;

    py::class_<QuinticPolyTrajectory>(m, "QuinticPolyTrajectory")
        .def(py::init<TFloat, const Vector6&, const Vector6&>(), py::arg("T"), py::arg("x_coeffs"), py::arg("y_coeffs"))
        .def(py::init<TFloat, const Vector3&, const Vector3&, const Vector3&, const Vector3&, bool>(),
            py::arg("T"), py::arg("x0"), py::arg("xT"), py::arg("y0"), py::arg("yT"), py::arg("relax_sx") = false)
        .def("__str__", py::overload_cast<const QuinticPolyTrajectory&>(&to_string))
        .def("__repr__", [](const QuinticPolyTrajectory& t) {
            stringstream ss; ss << "<QuinticPolyTrajectory with T=" << t.T() << "s>"; return ss.str();
        })
        .def_property_readonly("x_coeffs", &QuinticPolyTrajectory::x_coeffs)
        .def_property_readonly("y_coeffs", &QuinticPolyTrajectory::y_coeffs)
        .def_property_readonly("T", &QuinticPolyTrajectory::T)
        .def("point_at", &QuinticPolyTrajectory::point_at)
        .def("tangent_at", &QuinticPolyTrajectory::tangent_at)
        .def("velocity_at", &QuinticPolyTrajectory::velocity_at)
        .def("acceleration_at", &QuinticPolyTrajectory::acceleration_at)
        .def("periodize", &QuinticPolyTrajectory::periodize)
        ;

    py::class_<CTRATrajectory>(m, "CTRATrajectory")
        .def(py::init<TFloat, const Vector6&>(), py::arg("T"), py::arg("init_state"))
        .def(py::init<TFloat, Point, TFloat, TFloat, TFloat, TFloat>(), py::arg("T"),
            py::arg("p"), py::arg("theta"), py::arg("v"), py::arg("a") = 0, py::arg("omega") = 0)
        .def_property_readonly("initial_state", &CTRATrajectory::initial_state)
        .def_property_readonly("T", &CTRATrajectory::T)
        .def("point_at", &CTRATrajectory::point_at)
        .def("tangent_at", &CTRATrajectory::tangent_at)
        .def("velocity_at", &CTRATrajectory::velocity_at)
        // .def("acceleration_at", &CTRATrajectory::acceleration_at)
        .def("periodize", &CTRATrajectory::periodize)
        ;

    py::class_<HeteroSegment>(m, "HeteroSegment")
        .def(py::init<SegmentType, const vector<TFloat>&>())
        .def_readwrite("type", &HeteroSegment::type)
        .def_readwrite("params", &HeteroSegment::params)
        ;

    py::class_<HeteroPath>(m, "HeteroPath")
        .def(py::init<>())
        .def(py::init<const HeteroPath&>())
        .def(py::init<const Path&>())
        .def(py::init<const vector<Point>&, const vector<HeteroSegment>&, TFloat>(),
            py::arg("points"), py::arg("segments"), py::arg("s0") = 0)
        .def("__len__", &HeteroPath::size)
        .def_property_readonly("vertices", py::overload_cast<>(&HeteroPath::vertices, py::const_))
        .def_property_readonly("length", &HeteroPath::length)
        .def_property_readonly("segment_lengths", &HeteroPath::segment_lengths)
        .def("point_from", &HeteroPath::point_from)
        .def("point_at", &HeteroPath::point_at)
        .def("tangent_from", &HeteroPath::tangent_from)
        .def("tangent_at", &HeteroPath::tangent_at)
        .def("rasterize", &HeteroPath::rasterize, py::arg("resolution") = 0)
        // .def("densify", &HeteroPath::densify)
        ;

    py::module m_frenet = m.def_submodule("frenet", "Frenet related functionalitys");
    m_frenet.def("from_cartesian", py::overload_cast<const Path&, const Trajectory&>(frenet::from_cartesian));
    m_frenet.def("to_cartesian", py::overload_cast<const Path&, const Trajectory&>(frenet::to_cartesian));
    m_frenet.def("from_cartesian", py::overload_cast<const Path&, const Point&>(frenet::from_cartesian));
    m_frenet.def("to_cartesian", py::overload_cast<const Path&, const Point&>(frenet::to_cartesian));
    m_frenet.def("from_cartesian", py::overload_cast<const Path&, const Path&>(frenet::from_cartesian));
    m_frenet.def("to_cartesian", py::overload_cast<const Path&, const Path&>(frenet::to_cartesian));
    m_frenet.def("to_cartesian_rescale", &frenet::to_cartesian_rescale);

    m.def("distance", py::overload_cast<const Point&, const Point&>(traji::distance));
    m.def("distance", py::overload_cast<const Path&, const Point&>(traji::distance));
    m.def("arg_distance", py::overload_cast<const Path&, const Point&>(traji::arg_distance));
    m.def("tdistance", py::overload_cast<const Trajectory&, const Trajectory&>(traji::tdistance));
    m.def("intersection", py::overload_cast<const Path&, const Path&>(traji::intersection));
    m.def("arg_intersection", py::overload_cast<const Path&, const Path&>(traji::arg_intersection));
}
