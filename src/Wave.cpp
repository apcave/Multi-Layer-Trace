#include "Wave.hpp"

float Wave::omega = 0.0f;

bool Wave::operator==( Wave& other)  {
    if (type != other.type) {
        return false; // Different types cannot be equal
    }
    if (abs(angle - other.angle) >= 1e-6f) {
        return false; // Angles are not equal within a small tolerance
    }
    return true;
}

void Wave::clear() {
    p = std::complex<float>(0.0f,0.0f);
}

void Wave::accumulate( Wave& other) {
    p += other.p;
}

void Wave::print() const {
    if (abs(p) == 0.0f) {
        return; // Do not print if pressure is zero
    }
    std::cout << "Wave Type: " << (type == Type::P ? "P" : "S") 
                << ", Pressure: " << p 
                << ", Angle: " << angle
                << ", theta: " << angle.real() * 180.0f / M_PI << " degrees"
                << ", mag: " << std::abs(p)
                << std::endl;

}

bool Wave::isSurfaceWave() {
    if (abs(angle.real() - (M_PI / 2.0f)) < 1e-6f &&
        (abs(angle.imag()) < 1e-6f)) {
        std::cout << "Surface wave detected with angle: " << angle.real() * 180.0f / M_PI << " degrees" << std::endl;
        print();
        return true; // Surface wave if angle is close to ±90 degrees
    }
    return false; // Not a surface wave
}