#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>

#include "Composite.hpp"
#include "Medium.hpp"
#include "Wave.hpp"


namespace py = pybind11;

Composite comp;

// Note Angles can be represented as complex number and evanesce is supported.

std::vector<std::tuple<float, float, float, float>> getRp() {
    auto waves = comp.getRp();
    std::vector<std::tuple<float, float, float, float>> result;
    for (const auto& wave : waves) {
        result.emplace_back(wave.p.real(),wave.p.imag(), wave.angle.real(), wave.angle.imag());
    }
    return result;
}

std::vector<std::tuple<float, float, float, float>> getRs() {
    auto waves = comp.getRs();
    std::vector<std::tuple<float, float, float, float>> result;
    for (const auto& wave : waves) {
        result.emplace_back(wave.p.real(),wave.p.imag(), wave.angle.real(), wave.angle.imag());
    }
    return result;
}

std::vector<std::tuple<float, float, float, float>> getTp() {
    auto waves = comp.getTp();
    std::vector<std::tuple<float, float, float, float>> result;
    for (const auto& wave : waves) {
        result.emplace_back(wave.p.real(),wave.p.imag(), wave.angle.real(), wave.angle.imag());
    }
    return result;
}

std::vector<std::tuple<float, float, float, float>> getTs() {
    auto waves = comp.getTs();
    std::vector<std::tuple<float, float, float, float>> result;
    for (const auto& wave : waves) {
        result.emplace_back(wave.p.real(),wave.p.imag(), wave.angle.real(), wave.angle.imag());
    }
    return result;
}

void makeComposite( std::vector<float>& thickness,  std::vector<float>& density,
                    std::vector<float>& cp,  std::vector<float>& cs,
                    std::vector<float>& att_p,  std::vector<float>& att_s)
{
    std::cout << "Creating composite with " << thickness.size() << " layers." << std::endl;
    comp.makeComposite(thickness, density, cp, cs, att_p, att_s);
}

void properateWave( float p_real, float p_image, float angle, int type, float frequency)
{
    if (type < 0 || type > 1) {
        throw std::invalid_argument("Type must be 0 (P-wave) or 1 (S-wave).");
    }
    // Convert type to Wave::Type enum
    Wave::Type waveType = (type == 0) ? Wave::Type::P : Wave::Type::S;

    Wave::omega = frequency * 2.0f * M_PI; // Set the angular frequency

    std::complex<float> p = std::complex<float>(p_real, p_image); // Create complex pressure amplitude

    Wave incidentWave(waveType, p, angle);
    comp.properateWave(incidentWave);
}

PYBIND11_MODULE(composite_wave_response, m) {
    m.doc() = "Python interface for the Composite wave propagation system";

    m.def("makeComposite", &makeComposite, "Create the composite structure");
    m.def("properateWave", &properateWave, "Propagate a wave through the composite");
    m.def("getRp", &getRp, "Get reflected P-waves");
    m.def("getRs", &getRs, "Get reflected S-waves");
    m.def("getTp", &getTp, "Get transmitted P-waves");
    m.def("getTs", &getTs, "Get transmitted S-waves");
}