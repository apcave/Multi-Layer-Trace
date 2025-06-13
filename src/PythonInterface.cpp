#pragmo once
#include Composite.hpp

#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(composite_wave_response, m) {
    m.doc() = "Python interface for the Composite wave propagation system";

      py::class_<Wave>(m, "Wave")
        .def(py::init<>())
        .def_readwrite("amplitude", &Wave::amplitude)
        .def_readwrite("angle", &Wave::angle);

    pybind11::class_<Composite>(m, "Composite")
        .def(pybind11::init<>())
        .def("addMaterial", &Composite::addMaterial)
        .def("properateWave", &Composite::properateWave)
        .def("getRp", &Composite::getRp)
        .def("getRs", &Composite::getRs)
        .def("getTp", &Composite::getTp)
        .def("getTs", &Composite::getTs);
}