#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include "traji.hpp"

namespace py = pybind11;
using namespace traji;
using namespace std;

PYBIND11_MODULE(_lib2, m) {

    py::class_<Point>(m, "Point")
        .def(py::init<TFloat, TFloat>())
        .def(py::init<Point>())
        .def(py::init([](py::object point) {
            if (py::hasattr(point, "x") && py::hasattr(point, "y"))
                return Point(py::cast<TFloat>(point.attr("x")), py::cast<TFloat>(point.attr("y")));
            else
            {
                vector<TFloat> data = py::cast<vector<TFloat>>(point);
                if (data.size() != 2)
                    throw py::value_error("The length of point values should be 2!");
                return Point(data[0], data[1]);
            }
        }))
        .def_property_readonly("x", [](const Point &p) { return p.get<0>(); })
        .def_property_readonly("y", [](const Point &p) { return p.get<1>(); })
        .def("__str__", [](const Point &p) { return to_string(p); })
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
        .def(py::init<size_t, size_t, TFloat>(), py::arg("component"), py::arg("segment"), py::arg("fraction"))
        .def(py::init<size_t, TFloat>(), py::arg("segment"), py::arg("fraction"))
        .def("__str__", [](const PathPosition &pos) { return to_string(pos); })
        .def("__repr__", [](const PathPosition &pos) { 
            stringstream ss; ss << "<PathPosition ";
            ss << pos.component << ", ";
            ss << pos.segment << ", ";
            ss << pos.fraction << ">";
            return ss.str();
        })
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("to_s", static_cast<TFloat (PathPosition::*)(const Path&)>(&PathPosition::to_s))
        // .def("to_s", static_cast<TFloat (PathPosition::*)(const HeteroPath&)>(&PathPosition::to_s))
        .def_static("from_s", static_cast<PathPosition (*)(const Path&, TFloat)>(PathPosition::from_s))
        // .def_static("from_s", static_cast<PathPosition (*)(const HeteroPath&, TFloat)>(PathPosition::from_s))
        .def_static("from_s", static_cast<vector<PathPosition> (*)(const Path&, const vector<TFloat>&)>(PathPosition::from_s))
        .def("to_t", &PathPosition::to_t)
        .def_static("from_t", static_cast<PathPosition (*)(const Trajectory&, TFloat)>(PathPosition::from_t))
        .def_static("from_t", static_cast<vector<PathPosition> (*)(const Trajectory&, const vector<TFloat>&)>(PathPosition::from_t))
        ;

    py::class_<Path>(m, "Path")
        .def(py::init<>())
        .def(py::init<const vector<Point>&>())
        .def(py::init<const vector<Point>&, TFloat>())
        .def(py::init([](py::object point) {
        }
        ;
        
    py::class_<HeteroPath>(m, "HeteroPath")
        ;
}
