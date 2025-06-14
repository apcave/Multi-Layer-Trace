#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>

#include "Composite.hpp"
#include "Medium.hpp"
#include "Wave.hpp"


namespace py = pybind11;

Composite comp;

std::vector<Wave> getRp() {
    // Return an empty vector or your actual data
    return comp.getRp();
}

std::vector<Wave> getRs() {
    return comp.getRs();
}

std::vector<Wave> getTp() {
    return comp.getTp();
}

std::vector<Wave> getTs() {
    return comp.getTs();
}

void makeComposite( std::vector<float>& thickness,  std::vector<float>& density,
                    std::vector<float>& cp,  std::vector<float>& cs,
                    std::vector<float>& att_p,  std::vector<float>& att_s)
{
    comp.makeComposite(thickness, density, cp, cs, att_p, att_s);
}

void properateWave( Wave& incidentWave)
{
    comp.properateWave(incidentWave);
}

PYBIND11_MODULE(composite_wave_response, m) {
    m.doc() = "Python interface for the Composite wave propagation system";

      py::class_<Wave>(m, "Wave")
        .def(py::init<>())
        .def_readwrite("amplitude", &Wave::p)
        .def_readwrite("angle", &Wave::angle);

    pybind11::class_<Composite>(m, "Composite")
        .def(pybind11::init<>())
        .def("makeComposite", &Composite::makeComposite)
        .def("properateWave", &Composite::properateWave)
        .def("getRp", &Composite::getRp)
        .def("getRs", &Composite::getRs)
        .def("getTp", &Composite::getTp)
        .def("getTs", &Composite::getTs);
}